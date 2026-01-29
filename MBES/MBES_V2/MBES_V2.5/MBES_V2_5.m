%% 使用Field II计算多波束测深声纳(MBES) 
% 核心修正：
% 1. [物理] Mills Cross: 发射阵(Tx)沿航向(X)排列以产生跨航向扇面；接收阵(Rx)沿跨航向(Y)排列。
% 2. [逻辑] 单次Ping仿真: 模拟真实的单次宽波束发射，而非多次扫描。
% 3. [算法] Delay-and-Sum: 手动实现频域/时域波束形成。

clear all;
close all;

%% 1. 基础设置与常量
% 检查路径
if isempty(which('field_init'))
    addpath(genpath('./Field_II')); 
    addpath(genpath('./ustb')); 
end

c0 = 1500;      % 声速 [m/s]
fs = 4e6;       % 采样频率 [Hz] (降低采样率以节省内存)
dt = 1/fs;      

field_init(0);
set_field('c', c0);
set_field('fs', fs);
set_field('use_rectangles', 1);

%% 2. 换能器定义 (修正后的 Mills Cross)
% 坐标系: x(沿航向), y(跨航向), z(深度)

f0 = 200e3;           
lambda = c0/f0;       
kerf = 0.05*lambda;   

% === 发射阵列 (Tx) - 必须沿 X轴 (沿航向) 排列 ===
% 原理：X轴孔径长 -> X方向波束窄；Y轴孔径短 -> Y方向波束宽(扇面)
Tx_N = 64;                          
Tx_pitch = lambda/2;                
Tx_width = Tx_pitch - kerf; 
Tx_height = 15e-3; % Y方向高度

% Field II 默认 linear_array 沿 X 轴 -> 正确
Th = xdc_linear_array(Tx_N, Tx_width, Tx_height, kerf, 1, 1, [0 0 Inf]);

% === 接收阵列 (Rx) - 必须沿 Y轴 (跨航向) 排列 ===
% 原理：Y轴孔径长 -> 可在跨航向区分角度
Rx_N = 128;                          
Rx_pitch = lambda/2;                 
Rx_width = 15e-3;            % X方向宽度
Rx_height = Rx_pitch - kerf; % Y方向高度 (作为阵元长边)

% 手动定义沿 Y 轴的阵列
rect = zeros(Rx_N, 19);
center_rx = zeros(Rx_N, 3);
y_positions = ((1:Rx_N) - (Rx_N+1)/2) * Rx_pitch;

for i = 1:Rx_N
    cx = 0; cy = y_positions(i); cz = 0;
    rect(i, 1) = i;
    rect(i, 2:4)   = [-Rx_width/2, cy-Rx_height/2, 0];  
    rect(i, 5:7)   = [+Rx_width/2, cy-Rx_height/2, 0];  
    rect(i, 8:10)  = [+Rx_width/2, cy+Rx_height/2, 0];  
    rect(i, 11:13) = [-Rx_width/2, cy+Rx_height/2, 0];  
    rect(i, 14) = 1.0; 
    rect(i, 15) = Rx_width; 
    rect(i, 16) = Rx_height;   
    rect(i, 17:19) = [cx, cy, cz];
    center_rx(i, :) = [cx, cy, cz];
end
Rh = xdc_rectangles(rect, center_rx, [0 0 Inf]);

disp('========================================');
disp('MBES换能器配置 (Mills Cross - 修正版)');
disp('========================================');
fprintf('频率: %.1f kHz\n', f0/1e3);
fprintf('Tx (沿X轴): %d 阵元 -> 产生跨航向扇面\n', Tx_N);
fprintf('Rx (沿Y轴): %d 阵元 -> 接收跨航向回波\n', Rx_N);
disp('========================================');

%% 3. 脉冲定义
pulse_len = 1e-3; 
BW = 0.5 * f0;
t_p = 0:dt:pulse_len;
excitation = chirp(t_p, f0-BW/2, pulse_len, f0+BW/2) .* hamming(numel(t_p))';
excitation = excitation - mean(excitation);

impulse_response = sin(2*pi*f0*(0:dt:2/f0));
impulse_response = impulse_response .* hanning(length(impulse_response))';

xdc_excitation(Th, excitation);
xdc_impulse(Th, impulse_response);
xdc_impulse(Rh, impulse_response);

% Tx 聚焦：不偏转，聚焦在正下方远场，利用阵元本身特性形成扇面
xdc_center_focus(Th, [0 0 0]);
xdc_focus(Th, 0, [0 0 1000]); 
xdc_apodization(Th, 0, ones(1, Tx_N));

% Rx 聚焦：设置为动态聚焦全开(0延时)，以便获取原始通道数据
xdc_center_focus(Rh, [0 0 0]);
xdc_focus(Rh, 0, [0 0 1000]);
xdc_apodization(Rh, 0, ones(1, Rx_N));

%% 4. 海底地形生成
seabed_length_along = 40;     
seabed_width_across = 80;     
seabed_depth_mean = 50;       
N_seabed_scatter = 8000;       

Nx = round(sqrt(N_seabed_scatter * seabed_length_along / seabed_width_across));
Ny = round(N_seabed_scatter / Nx);
x_sb = linspace(-seabed_length_along/2, seabed_length_along/2, Nx);
y_sb = linspace(-seabed_width_across/2, seabed_width_across/2, Ny);
[Xg, Yg] = meshgrid(x_sb, y_sb);

% 地形函数
Zg = seabed_depth_mean + 1.5*sin(2*pi*Yg/30) + 1.0*cos(2*pi*Xg/20) + 0.3*randn(size(Xg));
% 添加小山丘
Zg = Zg - 12 * exp(-((Xg).^2 + (Yg-10).^2)/80); 

pos_seabed = [Xg(:), Yg(:), Zg(:)];
amp_seabed = 0.1 + rand(size(pos_seabed,1),1)*0.5;

%% 5. [参考 V2] 场景可视化
disp('正在生成场景预览...');
figure('Name', 'MBES场景预览', 'Position', [50 60 1200 600]);

% 子图1: 3D视图
subplot(1,2,1);
scatter3(pos_seabed(1:20:end,1), pos_seabed(1:20:end,2), pos_seabed(1:20:end,3), ...
         5, pos_seabed(1:20:end,3), 'filled'); 
hold on;
% 换能器
plot3(0, 0, 0, 'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'r');
% 绘制发射阵列方向 (Tx沿X轴 - 红色)
quiver3(0, 0, 0, 15, 0, 0, 'r', 'LineWidth', 3, 'MaxHeadSize', 0.5);
text(15,0,0, 'Tx (沿航向)', 'Color', 'r');
% 绘制接收阵列方向 (Rx沿Y轴 - 蓝色)
quiver3(0, 0, 0, 0, 15, 0, 'b', 'LineWidth', 3, 'MaxHeadSize', 0.5);
text(0,15,0, 'Rx (跨航向)', 'Color', 'b');

set(gca, 'ZDir', 'reverse');
axis equal tight; grid on;
xlabel('沿航向 x [m]'); ylabel('跨航向 y [m]'); zlabel('深度 z [m]');
title('3D场景视图 (修正版配置)');
view(35, 25); colormap(parula);

% 子图2: 跨航向覆盖示意
subplot(1,2,2);
swath_angle = 60;
range_tx = 70;
% 绘制波束扇面
theta = linspace(-swath_angle, swath_angle, 30);
y_fan = range_tx * sind(theta);
z_fan = range_tx * cosd(theta);
patch([0, y_fan, 0], [0, z_fan, 0], 'c', 'FaceAlpha', 0.1, 'EdgeColor', 'b');
hold on;
% 绘制海底剖面(x=0)
idx_prof = abs(pos_seabed(:,1)) < 2;
scatter(pos_seabed(idx_prof,2), pos_seabed(idx_prof,3), 10, 'k', 'filled');
set(gca, 'YDir', 'reverse');
axis equal tight; grid on;
xlabel('跨航向 y [m]'); ylabel('深度 z [m]');
title('跨航向波束覆盖示意');
drawnow;

%% 6. 仿真计算 (Single Ping)
disp(' ');
disp('========================================');
disp('开始Field II仿真 (Single Ping)...');
disp('========================================');

tic;
% Field II计算核心 - 一次性计算所有通道回波
[rf_data, t_start] = calc_scat_multi(Th, Rh, pos_seabed, amp_seabed);
calc_time = toc;
disp(['声场计算完成! 用时: ', num2str(calc_time, '%.1f'), ' 秒']);

%% 7. 信号处理与波束形成
disp('正在进行波束形成 (DAS)...');

% 匹配滤波
MF = flipud(excitation(:)); 
rf_filtered = zeros(size(rf_data));
for i = 1:Rx_N
    rf_filtered(:,i) = conv(rf_data(:,i), MF, 'same');
end

% 7.2 波束形成
angle_scan = linspace(-60, 60, 121); % 高分辨率扫描
Na = length(angle_scan);
depth_axis = linspace(seabed_depth_mean-20, seabed_depth_mean+20, 512);
n_depth = length(depth_axis);

bf_image = zeros(n_depth, Na);
bottom_detect_r = zeros(Na, 1);
bottom_detect_z = zeros(Na, 1);

time_axis = t_start + (0:size(rf_filtered,1)-1)' * dt;
rx_pos_y = center_rx(:, 2)'; % 接收阵元Y坐标

wb = waitbar(0, '波束形成中...');
for n = 1:Na
    waitbar(n/Na, wb);
    theta = angle_scan(n);
    
    % 计算延迟: delay = y * sin(theta) / c
    delay_dist = rx_pos_y * sind(theta);
    delay_time = delay_dist / c0;
    
    % 简单的延迟求和 (此处为了速度使用最近邻/线性插值近似)
    % 实际中可用频域相移
    beam_trace = zeros(size(time_axis));
    for ch = 1:Rx_N
        % 将 t 时刻的数据 搬移到 t - delay 时刻
        % 也就是 t_new = t_old - delay
        t_query = time_axis - delay_time(ch);
        sig_ch = interp1(time_axis, rf_filtered(:, ch), t_query, 'linear', 0);
        beam_trace = beam_trace + sig_ch;
    end
    
    env = abs(beam_trace);
    
    % TVG
    r_val = time_axis * c0 / 2;
    env = env .* (r_val .^ 1.5);
    
    % 转换到深度域显示
    z_val_beam = r_val * cosd(theta);
    bf_image(:, n) = interp1(z_val_beam, env, depth_axis, 'linear', 0);
    
    % 底检测
    [~, idx] = max(env(r_val > 20));
    if ~isempty(idx)
        valid_r = r_val(r_val > 20);
        r_det = valid_r(idx);
        bottom_detect_r(n) = r_det;
        bottom_detect_z(n) = r_det * cosd(theta);
    else
        bottom_detect_r(n) = NaN;
    end
end
close(wb);

% 归一化dB
bf_image_db = 20*log10(bf_image / max(bf_image(:)) + eps);
bf_image_db(bf_image_db < -60) = -60;

%% 8. [参考 V2] 结果显示与重建
disp('正在生成最终报告图...');

figure('Name', 'MBES仿真结果报告', 'Position', [100, 50, 1400, 800]);

% --- 子图1: 极坐标/扇形图像 ---
subplot(2,2,1);
imagesc(angle_scan, depth_axis, bf_image_db);
colormap(hot); colorbar;
clim([-50 0]);
xlabel('跨航向角度 [°]'); ylabel('深度 [m]');
title('波束形成图像 (Angle-Depth)');
set(gca, 'YDir', 'reverse');
grid on;

% --- 子图2: 直角坐标成像 ---
subplot(2,2,2);
[ANGLE_grid, DEPTH_grid] = meshgrid(angle_scan*pi/180, depth_axis);
Y_grid = DEPTH_grid .* sin(ANGLE_grid);
Z_grid = DEPTH_grid .* cos(ANGLE_grid);

pcolor(Y_grid, Z_grid, bf_image_db);
shading interp;
clim([-50 0]);
colorbar;
xlabel('跨航向距离 Y [m]'); ylabel('深度 Z [m]');
title('声纳图像 (笛卡尔坐标)');
set(gca, 'YDir', 'reverse');
axis equal tight;
grid on;

% --- 子图3: 3D点云重建 ---
subplot(2,2,3);
valid = bottom_detect_r > 0 & ~isnan(bottom_detect_r);
ang_valid = angle_scan(valid);
r_valid = bottom_detect_r(valid);

y_rec = r_valid(:) .* sind(ang_valid(:));
z_rec = r_valid(:) .* cosd(ang_valid(:));
x_rec = zeros(size(y_rec)); % 单ping假设x=0

scatter3(x_rec, y_rec, z_rec, 20, z_rec, 'filled');
hold on;
plot3(0,0,0, 'rp', 'MarkerSize',10, 'MarkerFaceColor','r');
set(gca, 'ZDir', 'reverse');
xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
title('3D 点云 (Single Ping)');
view(35, 25); grid on; axis equal;

% --- 子图4: 精度对比 ---
subplot(2,2,4);
% 真实剖面
idx_true = abs(pos_seabed(:,1)) < 1;
plot(pos_seabed(idx_true, 2), pos_seabed(idx_true, 3), 'k.', 'MarkerSize', 2);
hold on;
plot(y_rec, z_rec, 'b.-', 'LineWidth', 1.5);
set(gca, 'YDir', 'reverse');
legend('真实海底 (x=0)', '声纳检测', 'Location', 'SouthEast');
title('海底轮廓重建对比');
xlabel('跨航向 Y [m]'); ylabel('深度 Z [m]');
grid on; axis tight;

%% 9. 统计输出
disp(' ');
disp('========================================');
disp('仿真统计');
disp('========================================');
swath_w = max(y_rec) - min(y_rec);
fprintf('覆盖宽度: %.1f m (水深的 %.1f 倍)\n', swath_w, swath_w/seabed_depth_mean);
fprintf('有效波束: %d / %d\n', sum(valid), Na);
disp('========================================');

field_end;