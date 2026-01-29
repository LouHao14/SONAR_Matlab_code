%% MBES (多波束测深声纳) 完整仿真 - 物理修正版 (Corrected Mills Cross)
%
% 修正说明:
% 1. [物理架构] 交换 Tx/Rx 方向: Tx沿X轴(沿航向)，Rx沿Y轴(跨航向)。
%    这是真实MBES的标准配置，产生垂直于航向的扇面。
% 2. [逻辑修复] 降低深度门限至 20m，防止海底山丘(35m)被截断。
% 3. [可视化] 恢复完整的 3D 场景和阵列布局预览图。
% 4. [算法] 使用 Delay-and-Sum 波束形成逻辑。

clear all;
close all;

%% 1. 初始化
% 请确保 Field II 工具箱在路径中
addpath(genpath("./ustb"))     % 根据你的实际路径修改
addpath(genpath("./Field_II")) % 根据你的实际路径修改

field_init(0);

%% 2. 基础参数
c0 = 1500;      % 声速 [m/s]
fs = 4e6;       % 采样率 [Hz] (调试用4MHz，高精度可用10MHz)
dt = 1/fs;

set_field('c', c0);
set_field('fs', fs);
set_field('use_rectangles', 1);

%% 3. 换能器定义 (Mills Cross - 正确物理方向)
f0 = 200e3;           
lambda = c0/f0;       
kerf = 0.05*lambda;   

disp('配置换能器阵列...');
% === 发射阵列 (Tx) - 沿航向 (X轴) ===
% 目的: 在X方向形成窄波束(高沿航向分辨率)，在Y方向形成宽扇面(覆盖测宽)
Tx_N = 64;
Tx_pitch = lambda/2;
Tx_element_length_x = Tx_pitch - kerf; % 阵元在X方向的长度
Tx_element_width_y  = 15e-3;           % 阵元在Y方向的宽度 (很窄，保证Y方向波束宽)

% Field II 的 xdc_linear_array 默认就是沿 X 轴的
Th = xdc_linear_array(Tx_N, Tx_element_length_x, Tx_element_width_y, ...
                      kerf, 1, 1, [0 0 Inf]);

% === 接收阵列 (Rx) - 沿跨航向 (Y轴) ===
% 目的: 在Y方向有大孔径，以便进行跨航向波束形成(区分左右)
Rx_N = 128;
Rx_pitch = lambda/2;
Rx_element_length_y = Rx_pitch - kerf; % 阵元在Y方向的长度
Rx_element_width_x  = 15e-3;           % 阵元在X方向的宽度

% 手动构建沿 Y 轴的矩形阵列
rect = zeros(Rx_N, 19);
center_rx = zeros(Rx_N, 3);
y_positions = ((1:Rx_N) - (Rx_N+1)/2) * Rx_pitch;

for i = 1:Rx_N
    rect(i, 1) = i;
    cx = 0; cy = y_positions(i); cz = 0;
    
    % 定义矩形角点 (注意: 这是一个平铺在Z=0平面，长边沿Y轴的条状)
    hx = Rx_element_width_x / 2;  % X半宽
    hy = Rx_element_length_y / 2; % Y半长
    
    % 角点顺序: 逆时针
    rect(i, 2:4)   = [-hx, cy-hy, 0]; 
    rect(i, 5:7)   = [+hx, cy-hy, 0];
    rect(i, 8:10)  = [+hx, cy+hy, 0];
    rect(i, 11:13) = [-hx, cy+hy, 0];
    
    rect(i, 14) = 1.0; % Apodization
    rect(i, 15) = Rx_element_width_x;  % Width
    rect(i, 16) = Rx_element_length_y; % Height (Field II naming convention)
    rect(i, 17:19) = [cx, cy, cz];
    center_rx(i, :) = [cx, cy, cz];
end
Rh = xdc_rectangles(rect, center_rx, [0 0 Inf]);

% === 激励与脉冲响应 ===
fractional_bandwidth = 0.5;
BW = fractional_bandwidth * f0;
pulse_len = 1e-3;
t_p = 0:dt:pulse_len;
excitation = chirp(t_p, f0-BW/2, pulse_len, f0+BW/2) .* hamming(numel(t_p))';
excitation = excitation - mean(excitation);

impulse_response = 1; % 理想冲激响应
xdc_excitation(Th, excitation);
xdc_impulse(Th, impulse_response);
xdc_baffle(Th, 0); 
xdc_impulse(Rh, impulse_response);
xdc_baffle(Rh, 0); 

% === 初始聚焦设置 ===
% Tx: 发射平面波 (聚焦无穷远)，垂直向下
xdc_center_focus(Th, [0 0 0]);
xdc_focus(Th, 0, [0 0 1000]); 

% Rx: 接收动态聚焦 (后续在循环中设置或使用动态延迟)
xdc_center_focus(Rh, [0 0 0]);
xdc_focus(Rh, 0, [0 0 1000]);

%% 4. 海底地形生成 (带山丘)
disp('生成海底地形...');
seabed_mean = 50;
scan_width = 100;
N_scatter = 3000; % 散射点数量 (调试用3000，高精度用30000)

Nx = round(sqrt(N_scatter * 0.6)); % 纵横比调整
Ny = round(N_scatter / Nx);
x_s = linspace(-30, 30, Nx);
y_s = linspace(-scan_width/2, scan_width/2, Ny);
[Xg, Yg] = meshgrid(x_s, y_s);

% 地形函数
Zg = seabed_mean + 2*sin(2*pi*Yg/40) + 0.5*randn(size(Xg));
% 添加山丘 (高度15m => 顶端深度35m)
dist_hill = sqrt(Xg.^2 + Yg.^2);
Zg = Zg - 15 * exp(-(dist_hill/15).^2);

points = [Xg(:), Yg(:), Zg(:)];
amps = 0.1 + 0.2*rand(size(points,1), 1);

%% 5. [示意图绘制] 场景预览
disp('绘制场景预览...');
figure('Name', 'MBES仿真场景', 'Position', [50, 50, 1200, 700]);

% --- 3D 视图 ---
subplot(2,2,1);
sample_idx = 1:5:length(points); % 降采样绘图
scatter3(points(sample_idx,1), points(sample_idx,2), points(sample_idx,3), 5, points(sample_idx,3), 'filled');
hold on;
% 绘制阵列位置
plot3(0,0,0, 'rp', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
% 绘制 Tx 方向 (红色箭头, 沿X轴)
quiver3(0,0,0, 10,0,0, 'r', 'LineWidth', 3, 'MaxHeadSize', 0.5);
% 绘制 Rx 方向 (蓝色箭头, 沿Y轴)
quiver3(0,0,0, 0,10,0, 'b', 'LineWidth', 3, 'MaxHeadSize', 0.5);
set(gca, 'ZDir', 'reverse'); grid on; axis equal; view(45, 30);
xlabel('X (沿航向)'); ylabel('Y (跨航向)'); zlabel('深度');
title('3D 场景与阵列方向');
legend('海底', '船', 'Tx阵列', 'Rx阵列');

% --- 顶视图 ---
subplot(2,2,2);
scatter(points(sample_idx,1), points(sample_idx,2), 5, points(sample_idx,3), 'filled');
hold on;
plot(0,0, 'rp', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
grid on; axis equal; xlim([-40 40]); ylim([-60 60]);
xlabel('X (沿航向)'); ylabel('Y (跨航向)');
title('顶视图 (颜色代表深度)');

% --- 剖面图 (x=0) ---
subplot(2,2,3);
profile_idx = abs(points(:,1)) < 2;
scatter(points(profile_idx,2), points(profile_idx,3), 10, 'b', 'filled');
set(gca, 'YDir', 'reverse'); grid on; axis tight;
ylim([20 70]);
xlabel('Y (跨航向) [m]'); ylabel('深度 [m]');
title('跨航向剖面 (X \approx 0)');

% --- 阵列布局放大图 ---
subplot(2,2,4);
% 画 Tx (X轴)
rectangle('Position', [-Tx_N*Tx_pitch*100/2, -0.5, Tx_N*Tx_pitch*100, 1], 'FaceColor', 'r');
hold on;
% 画 Rx (Y轴)
rectangle('Position', [-0.5, -Rx_N*Rx_pitch*100/2, 1, Rx_N*Rx_pitch*100], 'FaceColor', 'b');
xlabel('cm'); ylabel('cm'); title('Mills Cross 布局 (红Tx, 蓝Rx)');
axis equal; xlim([-40 40]); ylim([-40 40]);
grid on;

drawnow;

%% 6. 仿真与波束形成
disp('开始Field II计算 (Sector Scan Mode)...');
% 定义角度
angles = linspace(-60, 60, 61); 
Na = length(angles);

% 深度轴
depth_axis = linspace(20, 80, 512)';
bf_image = zeros(length(depth_axis), Na);
bottom_detect = zeros(Na, 1);

% 匹配滤波器
MF = conj(flipud(excitation(:)));
one_way_ir = conv(1, excitation);
two_way_ir = conv(one_way_ir, MF);
lag = length(two_way_ir)/2;

wb = waitbar(0, '正在扫描...');

% 注意：为了模拟的准确性，这里我们让 Rx 在每个角度聚焦
% Tx 保持垂直向下发射宽波束 (对于 Mills Cross，Tx沿X轴，自然在Y方向是宽的)
xdc_focus(Th, 0, [0 0 1000]); 

for n = 1:Na
    waitbar(n/Na, wb);
    angle = angles(n);
    
    % Rx 聚焦于该角度 (动态聚焦)
    % 聚焦点坐标: y = R*sin(theta), z = R*cos(theta)
    xdc_focus(Rh, 0, [0, 1000*sind(angle), 1000*cosd(angle)]);
    
    % 计算声场
    [v, t] = calc_scat_multi(Th, Rh, points, amps);
    
    if ~isempty(v)
        % 匹配滤波 & 简单的波束求和
        sig_mf = zeros(size(v));
        win = hamming(Rx_N);
        for m = 1:Rx_N
            sig_mf(:,m) = conv(v(:,m), MF, 'same') * win(m);
        end
        beam_sig = sum(sig_mf, 2); % DAS
        
        % 时间轴 -> 距离轴
        t_axis = ((0:length(beam_sig)-1)*dt + t - lag*dt)';
        r_axis = t_axis * c0 / 2;
        
        % TVG 补偿
        tvg = abs(r_axis).^1.5;
        beam_sig_comp = abs(beam_sig) .* tvg;
        
        % 成像插值
        d_from_r = r_axis * cosd(angle);
        % 防止插值出NaN导致图像全白
        col_data = interp1(d_from_r, abs(beam_sig), depth_axis, 'linear', 1e-12);
        bf_image(:,n) = col_data;
        
        % 底检测 (关键修正: 门限设为20m)
        min_detect_depth = 20; 
        valid = r_axis > (min_detect_depth / cosd(angle));
        
        if any(valid)
            sig_search = beam_sig_comp(valid);
            r_search = r_axis(valid);
            [~, max_idx] = max(sig_search);
            bottom_detect(n) = r_search(max_idx);
        else
            bottom_detect(n) = NaN;
        end
    else
        bf_image(:,n) = 1e-12;
        bottom_detect(n) = NaN;
    end
end
close(wb);

%% 7. 结果可视化
disp('显示结果...');
% dB 转换与归一化
img_log = 20*log10(bf_image + 1e-12);
img_log = img_log - max(img_log(:));
clim_range = [-50 0];

figure('Name', '仿真结果', 'Position', [100, 100, 1200, 600]);

% [左图] 声纳图像
subplot(1,2,1);
imagesc(angles, depth_axis, img_log);
colormap(hot); colorbar; clim(clim_range);
set(gca, 'YDir', 'reverse');
xlabel('跨航向角度 (deg)'); ylabel('深度 (m)');
title('波束-深度图像 (B-Scan)');

% [右图] 重建对比
subplot(1,2,2);
valid_b = ~isnan(bottom_detect);
ang_v = angles(valid_b);
rng_v = bottom_detect(valid_b);

% 极坐标转直角坐标
y_rec = rng_v(:) .* sind(ang_v(:)); % 使用(:)确保列向量
z_rec = rng_v(:) .* cosd(ang_v(:));

% 绘制真实切片
idx_slice = abs(points(:,1)) < 2;
plot(points(idx_slice,2), points(idx_slice,3), 'k.', 'MarkerSize', 2);
hold on;
plot(y_rec, z_rec, 'b.-', 'LineWidth', 1.5);

set(gca, 'YDir', 'reverse'); grid on; axis equal tight;
ylim([20 80]);
xlabel('跨航向距离 Y (m)'); ylabel('深度 Z (m)');
title('海底地形重建对比');
legend('真实地形', '声纳重建');

field_end;