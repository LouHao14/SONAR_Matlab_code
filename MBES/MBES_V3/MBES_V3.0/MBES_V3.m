%% 使用Field II计算多波束测深声纳(MBES)数据集并进行波束形成
%
% Mills Cross正确配置：
% - 发射阵列: 沿航向(x轴)排列 → 跨航向波束宽(形成扇形覆盖)
% - 接收阵列: 跨航向(y轴)排列 → 跨航向波束窄(实现角度分辨)
%
% 坐标系: x(沿航向), y(跨航向), z(深度向下)

%% 清除工作空间
clear all;
close all;

%% 添加工具箱路径
addpath(genpath("./ustb"))
addpath(genpath("./Field_II"))

%% 基本常量
c0 = 1500;      % 声速 [m/s]
fs = 4e6;      % 采样频率 [Hz]
dt = 1/fs;      % 采样步长 [s]

%% Field II 初始化
field_init(0);
set_field('c', c0);
set_field('fs', fs);
set_field('use_rectangles', 1);

%% 换能器参数
f0 = 200e3;           % 中心频率 [Hz]
lambda = c0/f0;       % 波长 [m]
kerf = 0.05*lambda;   % 阵元间隙 [m]

% ============================================
% 发射阵列 - 沿航向(x轴)排列
% 阵列长、跨航向波束宽 → 形成跨航向扇形覆盖
% ============================================
Tx_N = 32;                          % 发射阵元数量
Tx_pitch = lambda/2;                % 阵元间距 [m]
Tx_element_width = Tx_pitch - kerf; % 阵元宽度(x方向) [m]
Tx_element_height = 15e-3;          % 阵元高度(y方向) [m]

% ============================================
% 接收阵列 - 跨航向(y轴)排列
% 阵列长、跨航向波束窄 → 实现跨航向角度分辨
% ============================================
Rx_N = 128;                          % 接收阵元数量
Rx_pitch = lambda/2;                 % 阵元间距 [m]
Rx_element_width = 15e-3;            % 阵元宽度(x方向) [m]
Rx_element_height = Rx_pitch - kerf; % 阵元高度(y方向) [m]

disp('========================================');
disp('MBES换能器配置 (Mills Cross)');
disp('========================================');
fprintf('中心频率: %.1f kHz, 波长: %.2f mm\n', f0/1e3, lambda*1e3);
fprintf('发射阵列(沿x轴): %d 阵元, 长度 %.2f cm\n', Tx_N, Tx_N*Tx_pitch*1e2);
fprintf('接收阵列(沿y轴): %d 阵元, 长度 %.2f cm\n', Rx_N, Rx_N*Rx_pitch*1e2);
disp('========================================');

%% 场景参数
seabed_depth_mean = 50;        % 平均水深 [m]
seabed_roughness = 3;          % 海底起伏 [m]
seabed_length_along = 60;      % 沿航向范围 [m]
seabed_width_across = 100;     % 跨航向范围 [m]
N_seabed_scatter = 3e3;        % 散射点数量

% 跨航向角度(波束形成方向)
across_track_angles = linspace(-60, 60, 61);  % 更多角度以获得更好分辨率
Na = length(across_track_angles);

disp(' ');
fprintf('平均水深: %.0f m, 角度范围: [%.0f°, %.0f°], %d个波束\n', ...
        seabed_depth_mean, across_track_angles(1), across_track_angles(end), Na);

%% 脉冲定义
fractional_bandwidth = 0.5;
BW = fractional_bandwidth * f0;
pulse_len = 1e-3;

t_p = 0:dt:pulse_len;
excitation = chirp(t_p, f0-BW/2, pulse_len, f0+BW/2) .* hamming(numel(t_p))';
excitation = excitation - mean(excitation);

impulse_response = 1;
MF = conj(flipud(excitation(:)));
two_way_ir = conv(conv(impulse_response, excitation), MF);
lag = length(two_way_ir)/2 + 1;

% 显示脉冲
figure('Name', 'MBES信号分析');
plot((0:(length(two_way_ir)-1))*dt - lag*dt, two_way_ir);
hold on; grid on; axis tight;
plot((0:(length(two_way_ir)-1))*dt - lag*dt, abs(hilbert(two_way_ir)), 'r', 'LineWidth', 1.5);
plot([0 0], ylim, 'g--', 'LineWidth', 1.5);
legend('双程脉冲', '包络', '零延迟');
title('双程冲激响应'); xlabel('时间 [s]'); ylabel('幅度');

%% 创建换能器阵列
% ============================================
% 发射阵列 - 沿x轴排列 (使用标准xdc_linear_array)
% Field II默认沿x轴排列，正好符合发射阵要求
% ============================================
noSubAz_Tx = max(1, round(Tx_element_width/(lambda/8)));
noSubEl_Tx = max(1, round(Tx_element_height/(lambda/8)));
Th = xdc_linear_array(Tx_N, Tx_element_width, Tx_element_height, ...
                      kerf, noSubAz_Tx, noSubEl_Tx, [0 0 Inf]);

% ============================================
% 接收阵列 - 沿y轴排列 (使用xdc_rectangles手动构建)
% ============================================
rect = zeros(Rx_N, 19);
center_rx = zeros(Rx_N, 3);

% 接收阵元沿y轴排列
y_positions = ((1:Rx_N) - (Rx_N+1)/2) * Rx_pitch;

for i = 1:Rx_N
    rect(i, 1) = i;  % 阵元编号
    
    cx = 0;
    cy = y_positions(i);
    cz = 0;
    
    % 阵元尺寸: width沿x, height沿y
    hw = Rx_element_width/2;   % x方向半宽
    hh = Rx_element_height/2;  % y方向半高
    
    % 四角点(逆时针)
    rect(i, 2:4)   = [-hw, cy-hh, 0];
    rect(i, 5:7)   = [+hw, cy-hh, 0];
    rect(i, 8:10)  = [+hw, cy+hh, 0];
    rect(i, 11:13) = [-hw, cy+hh, 0];
    
    rect(i, 14) = 1.0;                    % 加权
    rect(i, 15) = Rx_element_width;       % width
    rect(i, 16) = Rx_element_height;      % height
    rect(i, 17:19) = [cx, cy, cz];        % 中心
    
    center_rx(i, :) = [cx, cy, cz];
end

Rh = xdc_rectangles(rect, center_rx, [0 0 Inf]);

% 设置激励和响应
xdc_excitation(Th, excitation);
xdc_impulse(Th, impulse_response);
xdc_baffle(Th, 0);
xdc_center_focus(Th, [0 0 0]);

xdc_impulse(Rh, impulse_response);
xdc_baffle(Rh, 0);
xdc_center_focus(Rh, [0 0 0]);

disp('换能器阵列创建完成');

%% 生成海底地形
disp('正在生成海底地形...');

Nx = round(sqrt(N_seabed_scatter * seabed_length_along / seabed_width_across));
Ny = round(N_seabed_scatter / Nx);

x_seabed = linspace(-seabed_length_along/2, seabed_length_along/2, Nx);
y_seabed = linspace(-seabed_width_across/2, seabed_width_across/2, Ny);
[Xg, Yg] = meshgrid(x_seabed, y_seabed);

% 海底地形
Zg = seabed_depth_mean + ...
     seabed_roughness * (sin(2*pi*Xg/30) + ...
                         1.2*sin(2*pi*Yg/40) + ...
                         0.8*sin(2*pi*Xg/10).*cos(2*pi*Yg/15) + ...
                         0.5*randn(size(Xg)));

% 中心山丘
hill_height = 15;
hill_radius = 15;
Zg = Zg - hill_height * exp(-((Xg.^2 + Yg.^2)/hill_radius^2));

point_position = [Xg(:), Yg(:), Zg(:)];
point_amplitudes = 0.1 + 0.2*abs(randn(size(point_position,1), 1));

fprintf('海底散射点: %d\n', size(point_position, 1));

%% 场景可视化
figure('Name', 'MBES场景预览', 'Position', [50 50 1400 600]);

subplot(1,3,1);
scatter3(point_position(1:50:end,1), point_position(1:50:end,2), ...
         point_position(1:50:end,3), 3, point_position(1:50:end,3), 'filled');
hold on;
plot3(0, 0, 0, 'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'r');
% 发射阵方向(x轴)
quiver3(0, 0, 0, 15, 0, 0, 'b', 'LineWidth', 3, 'MaxHeadSize', 1);
text(16, 0, 0, 'Tx(x)', 'FontSize', 12, 'Color', 'b');
% 接收阵方向(y轴)
quiver3(0, 0, 0, 0, 15, 0, 'r', 'LineWidth', 3, 'MaxHeadSize', 1);
text(0, 16, 0, 'Rx(y)', 'FontSize', 12, 'Color', 'r');
set(gca, 'ZDir', 'reverse');
axis equal tight; grid on;
xlabel('沿航向 x [m]'); ylabel('跨航向 y [m]'); zlabel('深度 z [m]');
title('3D场景 (Tx沿x轴, Rx沿y轴)');
view(35, 25); colorbar;

subplot(1,3,2);
idx_profile = abs(Xg(:)) < 2;
scatter(point_position(idx_profile,2), point_position(idx_profile,3), ...
        5, point_position(idx_profile,3), 'filled');
hold on;
plot(0, 0, 'rp', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
set(gca, 'YDir', 'reverse');
axis equal tight; grid on;
xlabel('跨航向 y [m]'); ylabel('深度 z [m]');
title('跨航向剖面 (x≈0) - 期望重建结果');
colorbar;

subplot(1,3,3);
hold on; grid on; axis equal;
% 发射阵元(沿x轴)
x_tx = ((1:Tx_N) - (Tx_N+1)/2) * Tx_pitch;
for i = 1:Tx_N
    rectangle('Position', [(x_tx(i)-Tx_element_width/2)*100, -Tx_element_height/2*100, ...
              Tx_element_width*100, Tx_element_height*100], ...
              'FaceColor', [0.5 0.5 1], 'EdgeColor', 'b');
end
% 接收阵元(沿y轴)
for i = 1:Rx_N
    rectangle('Position', [-Rx_element_width/2*100, (y_positions(i)-Rx_element_height/2)*100, ...
              Rx_element_width*100, Rx_element_height*100], ...
              'FaceColor', [1 0.5 0.5], 'EdgeColor', 'r');
end
xlabel('x [cm]'); ylabel('y [cm]');
title('Mills Cross布局');
legend({'Tx(沿x)', 'Rx(沿y)'}, 'Location', 'best');
xlim([-15 15]); ylim([-50 50]);

drawnow;

%% Field II 仿真
disp(' ');
disp('========================================');
disp('开始Field II仿真...');
disp('========================================');

time_start = tic;
wb = waitbar(0, 'MBES仿真计算中...');

% 深度轴
depth_axis = linspace(seabed_depth_mean-25, seabed_depth_mean+15, 512)';
n_depth = length(depth_axis);

% 存储所有通道数据用于后续波束形成
all_channel_data = cell(1, 1);  % 单次发射
start_time = 0;

% ============================================
% 单次宽波束发射，接收阵列记录所有通道
% ============================================
fprintf('执行宽波束发射...\n');

% 发射设置：不聚焦(平面波)，形成跨航向宽扇形
xdc_apodization(Th, 0, ones(1, Tx_N));
xdc_focus(Th, 0, [0 0 1000]);  % 远场聚焦=平面波

% 接收设置
xdc_apodization(Rh, 0, ones(1, Rx_N));
xdc_focus_times(Rh, 0, zeros(1, Rx_N));

% 计算散射
[v, t] = calc_scat_multi(Th, Rh, point_position, point_amplitudes);
start_time = t;

if isempty(v)
    error('无回波数据!');
end

fprintf('接收数据: %d samples × %d channels\n', size(v,1), size(v,2));

% 匹配滤波
sig_mf = zeros(size(v));
for m = 1:Rx_N
    sig_mf(:,m) = conv(v(:,m), MF, 'same');
end

all_channel_data{1} = sig_mf;

waitbar(0.3, wb, '波束形成中...');

% ============================================
% 波束形成 - 对接收阵列进行相位延迟求和
% ============================================
disp('执行波束形成...');

bf_image = zeros(n_depth, Na);
bottom_depths = zeros(Na, 1);

% 时间轴
n_samples = size(sig_mf, 1);
time_axis = (0:n_samples-1) * dt + start_time - lag*dt;

% 接收阵元y坐标
rx_y = y_positions(:);

for n = 1:Na
    waitbar(0.3 + 0.7*n/Na, wb);
    
    theta = across_track_angles(n) * pi/180;  % 弧度
    
    % ============================================
    % 相控阵波束形成
    % 对于跨航向角度theta，计算各阵元的相位延迟
    % 延迟 = y * sin(theta) / c
    % ============================================
    delays = rx_y * sin(theta) / c0;  % [Rx_N × 1]
    
    % 延迟后求和
    beam_signal = zeros(n_samples, 1);
    
    for m = 1:Rx_N
        delay_samples = round(delays(m) / dt);
        
        if delay_samples >= 0
            if delay_samples < n_samples
                beam_signal(1:end-delay_samples) = beam_signal(1:end-delay_samples) + ...
                    sig_mf(delay_samples+1:end, m);
            end
        else
            delay_samples = abs(delay_samples);
            if delay_samples < n_samples
                beam_signal(delay_samples+1:end) = beam_signal(delay_samples+1:end) + ...
                    sig_mf(1:end-delay_samples, m);
            end
        end
    end
    
    % 取包络
    beam_envelope = abs(hilbert(beam_signal));
    
    % 时间转斜距
    range_axis = time_axis * c0 / 2;
    
    % 斜距转垂直深度
    depth_from_range = range_axis * cos(theta);
    
    % 插值到统一深度轴
    valid_idx = depth_from_range > 0;
    if sum(valid_idx) > 10
        bf_image(:, n) = interp1(depth_from_range(valid_idx), ...
                                  beam_envelope(valid_idx), ...
                                  depth_axis, 'linear', 0);
    end
    
    % 海底检测(最大值)
    [~, max_idx] = max(beam_envelope);
    if max_idx > 0 && max_idx <= length(range_axis)
        bottom_depths(n) = range_axis(max_idx);
    end
    
    if mod(n, 10) == 0
        fprintf('  波束 %d/%d (%.0f°)\n', n, Na, across_track_angles(n));
    end
end

close(wb);
time_elapsed = toc(time_start);
fprintf('仿真完成! 用时: %.1f 秒\n', time_elapsed);

% 归一化
bf_image = bf_image / max(bf_image(:) + eps);

%% 显示结果
bf_image_db = 20*log10(bf_image + eps);
bf_image_db = bf_image_db - max(bf_image_db(:));
bf_image_db(bf_image_db < -60) = -60;

figure('Name', 'MBES成像结果', 'Position', [100 50 1400 600]);

subplot(1,2,1);
imagesc(across_track_angles, depth_axis, bf_image_db);
colormap(hot); colorbar; clim([-60 0]);
xlabel('跨航向角度 [°]'); ylabel('深度 [m]');
title('波束形成图像 (角度-深度)');
set(gca, 'YDir', 'reverse'); grid on;

subplot(1,2,2);
[ANGLE_grid, DEPTH_grid] = meshgrid(across_track_angles*pi/180, depth_axis);
Y_grid = DEPTH_grid .* sin(ANGLE_grid);
Z_grid = DEPTH_grid .* cos(ANGLE_grid);
pcolor(Y_grid, Z_grid, bf_image_db);
shading interp; colormap(hot); colorbar; clim([-60 0]);
xlabel('跨航向距离 [m]'); ylabel('深度 [m]');
title('成像结果 (直角坐标)');
axis equal tight; set(gca, 'YDir', 'reverse'); grid on;

%% 3D海底重建
figure('Name', 'MBES 3D海底重建', 'Position', [100 100 1400 700]);

subplot(2,2,1);
imagesc(across_track_angles, depth_axis, bf_image_db);
colormap(hot); colorbar; clim([-60 0]);
xlabel('跨航向角度 [°]'); ylabel('深度 [m]');
title('波束-深度图');
set(gca, 'YDir', 'reverse'); grid on;

subplot(2,2,2);
valid_idx = bottom_depths > 0;
plot(across_track_angles(valid_idx), bottom_depths(valid_idx), 'b.-', 'LineWidth', 1.5);
hold on;
plot(across_track_angles, seabed_depth_mean*ones(size(across_track_angles)), 'r--', 'LineWidth', 2);
grid on;
xlabel('跨航向角度 [°]'); ylabel('检测深度(斜距) [m]');
title('海底深度检测');
legend('检测深度', '平均深度');
set(gca, 'YDir', 'reverse'); axis tight;

subplot(2,2,3);
angles_valid = across_track_angles(valid_idx) * pi/180;
depths_valid = bottom_depths(valid_idx);
y_bottom = depths_valid(:) .* sin(angles_valid(:));
z_bottom = depths_valid(:) .* cos(angles_valid(:));

scatter3(zeros(size(y_bottom)), y_bottom, z_bottom, 30, z_bottom, 'filled');
hold on;
plot3(0, 0, 0, 'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'r');
xlabel('沿航向 x [m]'); ylabel('跨航向 y [m]'); zlabel('深度 z [m]');
title('3D测深点云');
set(gca, 'ZDir', 'reverse'); axis equal; grid on; view(35, 25);
colorbar;

subplot(2,2,4);
% 真实海底剖面
idx_true = abs(Xg(:)) < 2;
plot(point_position(idx_true, 2), point_position(idx_true, 3), 'k.', 'MarkerSize', 2);
hold on;
plot(y_bottom, z_bottom, 'b.-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('跨航向 y [m]'); ylabel('深度 z [m]');
title('海底重建对比');
legend('真实海底', '重建结果', 'Location', 'best');
set(gca, 'YDir', 'reverse'); grid on; axis tight;

%% 统计
disp(' ');
disp('========================================');
disp('成像质量评估');
disp('========================================');
swath_width = max(y_bottom) - min(y_bottom);
fprintf('覆盖宽度: %.1f m\n', swath_width);
fprintf('覆盖比: %.2f (宽度/水深)\n', swath_width/seabed_depth_mean);
fprintf('有效波束: %d/%d (%.1f%%)\n', sum(valid_idx), Na, 100*sum(valid_idx)/Na);
disp('========================================');

%% 清理
field_end;
disp('MBES仿真完成!');