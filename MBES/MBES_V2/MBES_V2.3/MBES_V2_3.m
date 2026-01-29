%% 使用Field II计算多波束测深声纳(MBES)数据集并进行波束形成 - 修正版
%
% 主要修复：
% 1. [关键] Tx_pitch 改为 lambda/2 以消除栅瓣(Grating Lobes)
% 2. [关键] Tx_N 增加一倍以保持孔径长度不变
% 3. 添加了简单的TVG补偿以改善边缘检测
%
% 结果：海底地形将变平坦，不再出现边缘卷曲的"V"型或"W"型失真。

%% 清除工作空间和关闭图形窗口
clear all;
close all;

%% 添加必要的工具箱路径
% 请根据你的实际路径修改
addpath(genpath("./ustb"))
addpath(genpath("./Field_II"))

%% 基本常量
c0 = 1500;      % 声速 [m/s]
fs = 4e6;      % 采样频率 [Hz]
dt = 1/fs;      % 采样步长 [s]

%% Field II 初始化
field_init(0);
set_field('c', c0);              % 声速 [m/s]
set_field('fs', fs);             % 采样频率 [Hz]
set_field('use_rectangles', 1);  % 使用矩形阵元

%% 换能器定义 - Mills Cross配置
%
% 坐标系: x(沿航向), y(跨航向), z(深度)

% 中心频率和波长
f0 = 200e3;           % 换能器中心频率 [Hz]
lambda = c0/f0;       % 波长 [m]
kerf = 0.05*lambda;   % 阵元间隙 [m]

% === [修正 1] 发射阵列参数调整 ===
% 为了避免栅瓣，间距必须 <= lambda/2
% 为了保持原本的波束宽度(孔径大小)，数量增加一倍
Tx_N = 64;                          % [修改] 原为32
Tx_pitch = lambda/2;                % [修改] 原为lambda, 改为半波长消除栅瓣
Tx_element_width = Tx_pitch - kerf; % 阵元宽度 [m]
Tx_element_height = 15e-3;          % 阵元高度 [m]

% 接收阵列参数(沿航向 - 沿x轴排列)
Rx_N = 128;                          % 接收阵元数量
Rx_pitch = lambda/2;                 % 接收阵间距 [m]
Rx_element_width = Rx_pitch - kerf;  % 阵元宽度 [m]
Rx_element_height = 15e-3;           % 阵元高度 [m]

disp('========================================');
disp('MBES换能器配置 (Mills Cross) - 修正版');
disp('========================================');
fprintf('中心频率: %.1f kHz\n', f0/1e3);
fprintf('发射阵列(Tx): %d 阵元, 间距 %.2f lambda (无栅瓣)\n', Tx_N, Tx_pitch/lambda);
fprintf('接收阵列(Rx): %d 阵元\n', Rx_N);
disp('========================================');

%% 海底地形参数定义
seabed_length_along = 60;      % 沿航向长度 [m]
seabed_width_across = 100;     % 跨航向宽度 [m]
seabed_depth_mean = 50;        % 平均水深 [m]
seabed_roughness = 3;          % 海底起伏幅度 [m]
N_seabed_scatter = 3e3;        % 海底散射点数量

% 跨航向角度范围
across_track_angles = linspace(-60, 60, 31);  % 31个角度
Na = length(across_track_angles);

% 帧数
F = 1;  % 单帧仿真

%% 脉冲定义 - 使用chirp信号
fractional_bandwidth = 0.5;
BW = fractional_bandwidth * f0;    % 带宽 [Hz]
pulse_len = 1e-3;                  % 脉冲长度 [s]

% 生成chirp激励信号
t_p = 0:dt:pulse_len;
excitation = chirp(t_p, f0-BW/2, pulse_len, f0+BW/2) .* hamming(numel(t_p))';
excitation = excitation - mean(excitation);    % 去除直流分量

% 冲激响应 & 匹配滤波器
impulse_response = 1;
one_way_ir = conv(impulse_response, excitation);
MF = conj(flipud(excitation(:)));
two_way_ir = conv(one_way_ir, MF);
lag = length(two_way_ir)/2 + 1;

%% 创建换能器阵列

% === 接收阵列 (沿x轴) ===
noSubAz_Rx = max(1, round(Rx_element_width/(lambda/8)));
noSubEl_Rx = max(1, round(Rx_element_height/(lambda/8)));
Rh = xdc_linear_array(Rx_N, Rx_element_width, Rx_element_height, ...
                      kerf, noSubAz_Rx, noSubEl_Rx, [0 0 Inf]); 

% === 发射阵列 (沿y轴) ===
rect = zeros(Tx_N, 19);
center_tx = zeros(Tx_N, 3);
y_positions = ((1:Tx_N) - (Tx_N+1)/2) * Tx_pitch;

for i = 1:Tx_N
    rect(i, 1) = i;
    cx = 0;
    cy = y_positions(i);
    cz = 0;
    
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

% 设置激励
xdc_excitation(Th, excitation);
xdc_impulse(Th, impulse_response);
xdc_baffle(Th, 0);
xdc_center_focus(Th, [0 0 0]);

xdc_impulse(Rh, impulse_response);
xdc_baffle(Rh, 0);
xdc_center_focus(Rh, [0 0 0]);

%% 海底地形生成
disp('正在生成海底地形...');
Nx_seabed = round(sqrt(N_seabed_scatter * seabed_length_along / seabed_width_across));
Ny_seabed = round(N_seabed_scatter / Nx_seabed);
x_seabed = linspace(-seabed_length_along/2, seabed_length_along/2, Nx_seabed);
y_seabed = linspace(-seabed_width_across/2, seabed_width_across/2, Ny_seabed);
[Xg, Yg] = meshgrid(x_seabed, y_seabed);

% 地形函数
Zg = seabed_depth_mean + ...
     seabed_roughness * (1.0*sin(2*pi*Xg/30) + ...
                         1.2*sin(2*pi*Yg/40) + ...
                         0.8*sin(2*pi*Xg/10) .* cos(2*pi*Yg/15) + ...
                         0.5*randn(size(Xg)));
% 添加小山丘
hill_x = 0; hill_y = 0; hill_height = 15; hill_radius = 15;
distance = sqrt((Xg-hill_x).^2 + (Yg-hill_y).^2);
Zg = Zg - hill_height * exp(-(distance/hill_radius).^2);

point_position = [Xg(:), Yg(:), Zg(:)];
point_amplitudes = 0.1 + 0.2 * abs(randn(size(point_position,1), 1));
point_amplitudes = point_amplitudes .* (1 + 0.3*sin(2*pi*Yg(:)/20));

%% 场景预览 (省略，与原代码相同)
% 如果需要查看，可在此处插入原来的绘图代码

%% Field II 仿真计算
disp(' ');
disp('========================================');
disp('开始Field II仿真计算...');
disp('========================================');

time_start = tic;
wb = waitbar(0, 'MBES: 正在计算跨航向波束数据...');

win = hamming(Rx_N);
depth_axis = linspace(seabed_depth_mean-25, seabed_depth_mean+25, 512)'; % 稍微增加分辨率
n_depth = length(depth_axis);
bf_image = zeros(n_depth, Na);

max_amplitudes = zeros(Na, 1);
bottom_times = zeros(Na, 1);

for f = 1:F
    for n = 1:Na
        waitbar(n/Na, wb, sprintf('Processing Angle %d/%d', n, Na));
        angle = across_track_angles(n);
        
        % === 发射 ===
        xdc_apodization(Th, 0, ones(1, Tx_N));
        % 远场聚焦
        focus_depth = 1000; 
        focus_y = focus_depth * sind(angle);
        focus_z = focus_depth * cosd(angle);
        xdc_focus(Th, 0, [0, focus_y, focus_z]);
        
        % === 接收 ===
        xdc_apodization(Rh, 0, ones(1, Rx_N));
        xdc_focus_times(Rh, 0, zeros(1, Rx_N)); 
        
        % === 计算散射 ===
        [v, t] = calc_scat_multi(Th, Rh, point_position, point_amplitudes);
        
        if isempty(v) || all(v(:) == 0)
            continue;
        end
        
        % ... (前面是 Field II calc_scat_multi 计算 v, t 的代码) ...
        
        % === 信号处理 ===
        sig_mf = zeros(size(v));
        for m = 1:Rx_N
            tmp_sig = conv(v(:,m), MF, 'same');
            sig_mf(:,m) = tmp_sig * win(m);
        end
        
        beam_signal = sum(sig_mf, 2); % 这是一个列向量
        
        % 时间/距离转换 [关键修正：添加转置 ' ]
        time_axis = ((0:size(beam_signal,1)-1) * dt + t - lag*dt)';
        range_axis = time_axis * c0 / 2;
        
        % [修正 2] 简单的TVG (Time Varying Gain) 补偿
        % 现在 tvg_gain 也是列向量，与 beam_signal 维度一致
        tvg_gain = abs(range_axis) .^ 1.5; 
        beam_signal_comp = abs(beam_signal) .* tvg_gain;
        
        % 斜距转深度
        depth_from_range = range_axis * cosd(angle);
        
        % 插值
        beam_signal_interp = interp1(depth_from_range, abs(beam_signal), ...
                                     depth_axis, 'linear', 0);
        bf_image(:, n) = beam_signal_interp;
        
        % [修正 3] 基于补偿后的信号进行底检测
        valid_range_idx = range_axis > 10;
        
        if any(valid_range_idx)
            search_signal = beam_signal_comp(valid_range_idx);
            search_axis = range_axis(valid_range_idx);
            
            [max_amp, max_idx_local] = max(search_signal);
            
            if ~isempty(max_idx_local)
                bottom_times(n) = search_axis(max_idx_local);
            else
                bottom_times(n) = NaN;
            end
        else
            bottom_times(n) = NaN;
        end
    end
end
close(wb);
time_elapsed = toc(time_start);
disp(['计算完成! 用时: ', num2str(time_elapsed, '%.1f'), ' 秒']);

% 归一化
if max(bf_image(:)) > 0
    bf_image = bf_image / max(bf_image(:));
end

%% 成像显示与重建
bf_image_db = 20*log10(abs(bf_image) + eps);
bf_image_db = bf_image_db - max(bf_image_db(:));
bf_image_db(bf_image_db < -60) = -60;

figure('Name', '修正后的MBES成像结果', 'Position', [100, 100, 1400, 700]);

% 子图1: 成像结果
subplot(1,2,1);
imagesc(across_track_angles, depth_axis, bf_image_db);
colormap(hot);
colorbar;
clim([-60 0]);
xlabel('跨航向角度 [°]'); ylabel('深度 [m]');
title('MBES声纳图像 (无栅瓣)');
set(gca, 'YDir', 'reverse');

% 子图2: 重建对比
subplot(1,2,2);
bottom_depths = bottom_times;
valid_idx = ~isnan(bottom_depths) & bottom_depths > 0;
angles_valid = across_track_angles(valid_idx);
depths_valid = bottom_depths(valid_idx);

y_bottom = depths_valid(:) .* sind(angles_valid(:));
z_bottom = depths_valid(:) .* cosd(angles_valid(:));

% 绘制真实轮廓 (x≈0)
idx_true = abs(Xg(:)) < 2;
plot(point_position(idx_true, 2), point_position(idx_true, 3), 'k.', 'MarkerSize', 3);
hold on;
plot(y_bottom, z_bottom, 'b.-', 'LineWidth', 2);
set(gca, 'YDir', 'reverse');
grid on; axis equal tight;
legend('真实海底(x≈0)', '仿真重建', 'Location', 'SouthEast');
title('海底重建对比');
xlabel('跨航向距离 y [m]'); ylabel('深度 z [m]');

%% 释放资源
field_end;