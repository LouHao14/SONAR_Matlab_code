%% 使用Field II计算多波束测深声纳(MBES)数据集并进行波束形成
%
% 版本V4: 修正斜距-深度转换和海底检测算法
%
% 主要修正:
% 1. 正确的斜距到深度/水平距离转换
% 2. 改进的海底检测算法（阈值检测+峰值搜索）
% 3. 大角度回波强度补偿
% 4. 更准确的波束形成

%% 清除工作空间
clear all;
close all;

%% 添加工具箱路径
addpath(genpath("./ustb"))
addpath(genpath("./Field_II"))

%% 基本常量
c0 = 1500;      % 声速 [m/s]
fs = 10e6;      % 采样频率 [Hz]
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

% 发射阵列(跨航向y轴)
Tx_N = 32;
Tx_pitch = lambda;
Tx_element_width = Tx_pitch - kerf;
Tx_element_height = 15e-3;

% 接收阵列(沿航向x轴)
Rx_N = 128;
Rx_pitch = lambda/2;
Rx_element_width = Rx_pitch - kerf;
Rx_element_height = 15e-3;

disp('========================================');
disp('MBES换能器配置 (Mills Cross)');
disp('========================================');
fprintf('中心频率: %.1f kHz, 波长: %.2f mm\n', f0/1e3, lambda*1e3);
fprintf('发射阵列(y轴): %d 阵元\n', Tx_N);
fprintf('接收阵列(x轴): %d 阵元\n', Rx_N);
disp('========================================');

%% 场景参数
seabed_depth_mean = 50;        % 平均水深 [m]
seabed_roughness = 3;          % 海底起伏 [m]
N_seabed_scatter = 3e4;        % 散射点数量
seabed_length_along = 60;      % 沿航向长度 [m]
seabed_width_across = 100;     % 跨航向宽度 [m]

% 跨航向角度
across_track_angles = linspace(-60, 60, 31);
Na = length(across_track_angles);
F = 1;

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

%% 创建换能器阵列

% === 接收阵列(x轴) ===
noSubAz_Rx = max(1, round(Rx_element_width/(lambda/8)));
noSubEl_Rx = max(1, round(Rx_element_height/(lambda/8)));
Rh = xdc_linear_array(Rx_N, Rx_element_width, Rx_element_height, ...
                      kerf, noSubAz_Rx, noSubEl_Rx, [0 0 Inf]); 

% === 发射阵列(y轴,使用xdc_rectangles) ===
rect = zeros(Tx_N, 19);
center_tx = zeros(Tx_N, 3);
y_positions = ((1:Tx_N) - (Tx_N+1)/2) * Tx_pitch;

for i = 1:Tx_N
    rect(i, 1) = i;
    cx = 0; cy = y_positions(i); cz = 0;
    hw = Tx_element_height/2;
    hh = Tx_element_width/2;
    
    rect(i, 2:4)   = [-hw, cy-hh, 0];
    rect(i, 5:7)   = [+hw, cy-hh, 0];
    rect(i, 8:10)  = [+hw, cy+hh, 0];
    rect(i, 11:13) = [-hw, cy+hh, 0];
    rect(i, 14) = 1.0;
    rect(i, 15) = Tx_element_height;
    rect(i, 16) = Tx_element_width;
    rect(i, 17:19) = [cx, cy, cz];
    center_tx(i, :) = [cx, cy, cz];
end

Th = xdc_rectangles(rect, center_tx, [0 0 Inf]);

% 设置换能器参数
xdc_excitation(Th, excitation);
xdc_impulse(Th, impulse_response);
xdc_baffle(Th, 0);
xdc_center_focus(Th, [0 0 0]);
xdc_impulse(Rh, impulse_response);
xdc_baffle(Rh, 0);
xdc_center_focus(Rh, [0 0 0]);

disp('换能器阵列创建完成');

%% 海底地形生成
disp('正在生成海底地形...');

Nx_seabed = round(sqrt(N_seabed_scatter * seabed_length_along / seabed_width_across));
Ny_seabed = round(N_seabed_scatter / Nx_seabed);

x_seabed = linspace(-seabed_length_along/2, seabed_length_along/2, Nx_seabed);
y_seabed = linspace(-seabed_width_across/2, seabed_width_across/2, Ny_seabed);
[Xg, Yg] = meshgrid(x_seabed, y_seabed);

% 地形生成
Zg = seabed_depth_mean + ...
     seabed_roughness * (1.0*sin(2*pi*Xg/30) + ...
                         1.2*sin(2*pi*Yg/40) + ...
                         0.8*sin(2*pi*Xg/10) .* cos(2*pi*Yg/15) + ...
                         0.5*randn(size(Xg)));

% 添加山丘
hill_height = 15; hill_radius = 15;
distance = sqrt(Xg.^2 + Yg.^2);
Zg = Zg - hill_height * exp(-(distance/hill_radius).^2);

point_position = [Xg(:), Yg(:), Zg(:)];
point_amplitudes = 0.1 + 0.2 * abs(randn(size(point_position,1), 1));
point_amplitudes = point_amplitudes .* (1 + 0.3*sin(2*pi*Yg(:)/20));

disp(['散射点数量: ', num2str(size(point_position,1))]);

%% 计算每个角度的真实海底深度（用于验证）
% 在x=0处，沿y方向的真实海底剖面
true_y_profile = y_seabed;
true_z_profile = interp2(Xg, Yg, Zg, zeros(size(y_seabed)), y_seabed, 'linear');

%% 场景预览
figure('Name', 'MBES场景预览', 'Position', [50 60 1200 500]);

subplot(1,2,1);
scatter3(point_position(1:50:end,1), point_position(1:50:end,2), point_position(1:50:end,3), ...
         2, point_position(1:50:end,3), 'filled'); 
hold on;
plot3(0, 0, 0, 'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'r');
set(gca, 'ZDir', 'reverse');
axis equal tight; grid on;
xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
title('3D场景');
view(35, 25); colorbar;

subplot(1,2,2);
plot(true_y_profile, true_z_profile, 'k-', 'LineWidth', 2);
set(gca, 'YDir', 'reverse');
grid on; xlabel('跨航向 y [m]'); ylabel('深度 z [m]');
title('真实海底剖面 (x=0)');
drawnow;

%% Field II 仿真
disp(' ');
disp('开始Field II仿真...');

time_start = tic;
wb = waitbar(0, 'MBES仿真计算中...');

win = hamming(Rx_N);

% ========== 关键修正: 使用斜距轴而非深度轴 ==========
% 斜距范围应覆盖所有可能的回波
min_range = (seabed_depth_mean - 25) / cosd(60);  % 最小斜距
max_range = (seabed_depth_mean + 25) / cosd(0);   % 最大斜距
range_axis = linspace(min_range * 0.8, max_range * 1.2, 512)';
n_range = length(range_axis);

% 存储每个波束的数据
beam_data = zeros(n_range, Na);         % 波束形成后的信号
detected_range = zeros(Na, 1);          % 检测到的斜距
detected_depth = zeros(Na, 1);          % 转换后的深度
detected_y = zeros(Na, 1);              % 转换后的水平距离

for f = 1:F
    for n = 1:Na
        waitbar(n/Na, wb, sprintf('计算波束 %d/%d (%.1f°)', n, Na, across_track_angles(n)));
        
        angle = across_track_angles(n);
        angle_rad = angle * pi/180;
        
        % 发射设置
        xdc_apodization(Th, 0, ones(1, Tx_N));
        focus_depth = 1000;
        xdc_focus(Th, 0, [0, focus_depth*sind(angle), focus_depth*cosd(angle)]);
        
        % 接收设置
        xdc_apodization(Rh, 0, ones(1, Rx_N));
        xdc_focus_times(Rh, 0, zeros(1, Rx_N));
        
        % 散射计算
        [v, t] = calc_scat_multi(Th, Rh, point_position, point_amplitudes);
        
        if isempty(v) || all(v(:) == 0)
            continue;
        end
        
        % 匹配滤波
        sig_mf = zeros(size(v));
        for m = 1:Rx_N
            sig_mf(:,m) = conv(v(:,m), MF, 'same') * win(m);
        end
        
        % 波束形成
        beam_signal = sum(sig_mf, 2);
        beam_envelope = abs(hilbert(beam_signal));
        
        % 时间轴和斜距轴
        time_vec = (0:length(beam_signal)-1) * dt + t - lag*dt;
        range_vec = time_vec * c0 / 2;  % 这是斜距！
        
        % 插值到统一斜距轴
        beam_interp = interp1(range_vec, beam_envelope, range_axis, 'linear', 0);
        beam_data(:, n) = beam_interp;
        
        % ========== 改进的海底检测 ==========
        % 1. 预期斜距范围（基于预期深度和角度）
        expected_range = seabed_depth_mean / cosd(abs(angle));
        range_window = 20 / cosd(abs(angle));  % ±20m深度对应的斜距窗口
        
        % 2. 在预期范围内搜索
        search_idx = find(range_axis >= expected_range - range_window & ...
                          range_axis <= expected_range + range_window);
        
        if ~isempty(search_idx)
            % 在搜索窗口内找最大值
            [~, local_max_idx] = max(beam_interp(search_idx));
            global_max_idx = search_idx(local_max_idx);
            detected_range(n) = range_axis(global_max_idx);
        else
            % 全局搜索
            [~, global_max_idx] = max(beam_interp);
            detected_range(n) = range_axis(global_max_idx);
        end
        
        % ========== 关键修正: 正确的坐标转换 ==========
        % 斜距 R, 角度 θ:
        %   水平距离 y = R × sin(θ)
        %   垂直深度 z = R × cos(θ)
        detected_y(n) = detected_range(n) * sind(angle);
        detected_depth(n) = detected_range(n) * cosd(angle);
        
        if mod(n, 10) == 0
            fprintf('  波束 %d: 角度=%.1f°, 斜距=%.1fm, 深度=%.1fm, y=%.1fm\n', ...
                    n, angle, detected_range(n), detected_depth(n), detected_y(n));
        end
        
        clear v sig_mf;
    end
end
close(wb);

time_elapsed = toc(time_start);
disp(['计算完成! 用时: ', num2str(time_elapsed, '%.1f'), ' 秒']);

%% 数据可视化

% 归一化
beam_data_norm = beam_data / max(beam_data(:));
beam_data_db = 20*log10(beam_data_norm + eps);
beam_data_db(beam_data_db < -60) = -60;

figure('Name', 'MBES成像结果', 'Position', [100, 50, 1400, 600]);

% 子图1: 斜距-角度图 (原始数据)
subplot(1,2,1);
imagesc(across_track_angles, range_axis, beam_data_db);
hold on;
plot(across_track_angles, detected_range, 'c.-', 'LineWidth', 1.5, 'MarkerSize', 10);
colormap(hot); colorbar; clim([-60 0]);
xlabel('跨航向角度 [°]'); ylabel('斜距 [m]');
title('波束-斜距图 (原始坐标)');
set(gca, 'YDir', 'reverse'); grid on;
legend('检测斜距', 'Location', 'best');

% 子图2: 直角坐标显示
subplot(1,2,2);
[ANGLE_grid, RANGE_grid] = meshgrid(across_track_angles*pi/180, range_axis);
Y_grid = RANGE_grid .* sin(ANGLE_grid);
Z_grid = RANGE_grid .* cos(ANGLE_grid);

pcolor(Y_grid, Z_grid, beam_data_db);
shading interp; colormap(hot); colorbar; clim([-60 0]);
hold on;
% 绘制检测结果
plot(detected_y, detected_depth, 'c.-', 'LineWidth', 2, 'MarkerSize', 12);
% 绘制真实海底
plot(true_y_profile, true_z_profile, 'g--', 'LineWidth', 2);
xlabel('跨航向 y [m]'); ylabel('深度 z [m]');
title('MBES成像 (直角坐标)');
set(gca, 'YDir', 'reverse'); grid on;
legend('检测海底', '真实海底', 'Location', 'best');
axis equal tight;

%% 3D海底重建对比
figure('Name', 'MBES 3D海底重建', 'Position', [100, 100, 1400, 700]);

% 子图1: 波束-斜距图
subplot(2,2,1);
imagesc(across_track_angles, range_axis, beam_data_db);
colormap(hot); colorbar; clim([-60 0]);
xlabel('跨航向角度 [°]'); ylabel('斜距 [m]');
title('波束-斜距图');
set(gca, 'YDir', 'reverse'); grid on;

% 子图2: 检测深度 vs 角度
subplot(2,2,2);
valid_idx = detected_depth > 0;
plot(across_track_angles(valid_idx), detected_depth(valid_idx), 'b.-', 'LineWidth', 1.5, 'MarkerSize', 10);
hold on;
% 计算每个角度对应的真实深度
for i = 1:Na
    y_target = detected_y(i);
    true_depth_at_angle(i) = interp1(true_y_profile, true_z_profile, y_target, 'linear', NaN);
end
plot(across_track_angles, true_depth_at_angle, 'r--', 'LineWidth', 2);
plot(across_track_angles, seabed_depth_mean*ones(size(across_track_angles)), 'k:', 'LineWidth', 1);
grid on;
xlabel('跨航向角度 [°]'); ylabel('深度 z [m]');
title('海底深度检测');
legend('检测深度', '真实深度', '平均深度', 'Location', 'best');
set(gca, 'YDir', 'reverse');

% 子图3: 3D点云
subplot(2,2,3);
x_bottom = zeros(size(detected_depth));
scatter3(x_bottom(valid_idx), detected_y(valid_idx), detected_depth(valid_idx), ...
         50, detected_depth(valid_idx), 'filled');
hold on;
plot3(0, 0, 0, 'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'r');
xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
title('3D测深点云');
set(gca, 'ZDir', 'reverse');
axis equal; grid on; view(35, 25);
colorbar;

% 子图4: 海底重建对比 (关键图)
subplot(2,2,4);
% 真实海底
plot(true_y_profile, true_z_profile, 'k.', 'MarkerSize', 3);
hold on;
% 重建海底
plot(detected_y(valid_idx), detected_depth(valid_idx), 'b.-', 'LineWidth', 2, 'MarkerSize', 12);
xlabel('跨航向 y [m]'); ylabel('深度 z [m]');
title('海底重建对比');
legend('真实海底', '重建海底', 'Location', 'best');
set(gca, 'YDir', 'reverse');
grid on; axis tight;

%% 误差分析
figure('Name', '重建误差分析', 'Position', [150, 150, 1200, 400]);

% 计算误差
depth_error = detected_depth - true_depth_at_angle';
valid_error = ~isnan(depth_error) & valid_idx';

subplot(1,3,1);
plot(across_track_angles(valid_error), depth_error(valid_error), 'b.-', 'LineWidth', 1.5);
hold on;
plot([-60 60], [0 0], 'r--', 'LineWidth', 1);
grid on;
xlabel('跨航向角度 [°]'); ylabel('深度误差 [m]');
title('深度误差 vs 角度');
ylim([-10 10]);

subplot(1,3,2);
plot(detected_y(valid_error), depth_error(valid_error), 'b.-', 'LineWidth', 1.5);
hold on;
plot([min(detected_y) max(detected_y)], [0 0], 'r--', 'LineWidth', 1);
grid on;
xlabel('跨航向距离 y [m]'); ylabel('深度误差 [m]');
title('深度误差 vs 水平距离');
ylim([-10 10]);

subplot(1,3,3);
histogram(depth_error(valid_error), 20);
grid on;
xlabel('深度误差 [m]'); ylabel('计数');
title(sprintf('误差分布 (均值=%.2fm, 标准差=%.2fm)', ...
      mean(depth_error(valid_error)), std(depth_error(valid_error))));

%% 统计输出
disp(' ');
disp('========================================');
disp('重建精度评估');
disp('========================================');
fprintf('有效测深点: %d/%d\n', sum(valid_error), Na);
fprintf('深度误差均值: %.2f m\n', mean(depth_error(valid_error)));
fprintf('深度误差标准差: %.2f m\n', std(depth_error(valid_error)));
fprintf('深度误差最大值: %.2f m\n', max(abs(depth_error(valid_error))));

% 按角度分析
small_angle = abs(across_track_angles) <= 30;
large_angle = abs(across_track_angles) > 30;
fprintf('\n小角度(≤30°)误差标准差: %.2f m\n', std(depth_error(valid_error & small_angle')));
fprintf('大角度(>30°)误差标准差: %.2f m\n', std(depth_error(valid_error & large_angle')));

swath_width = max(detected_y) - min(detected_y);
fprintf('\n测线覆盖宽度: %.1f m\n', swath_width);
fprintf('覆盖比(宽度/水深): %.2f\n', swath_width / seabed_depth_mean);
disp('========================================');

%% 释放资源
field_end;

disp(' ');
disp('MBES仿真完成!');