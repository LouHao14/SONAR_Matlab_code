%% 使用Field II计算多波束测深声纳(MBES)数据集并进行波束形成
%
% 本示例展示如何使用Field II创建标准船载多波束测深声纳(MBES)仿真
% 
% 主要功能:
% - 使用Mills Cross正交双阵列配置(发射阵+接收阵)
% - 生成真实的海底地形散射场景
% - 使用chirp信号进行声场仿真
% - 跨航向多波束接收
% - 3D海底地形重建
%
% 修正内容:
% - 使用xdc_rectangles正确实现90度旋转的发射阵列
% - 修正聚焦坐标设置
% - 修正变量定义问题
% - 完善波束形成流程

%% 清除工作空间和关闭图形窗口
clear all;
close all;

%% 添加必要的工具箱路径
addpath(genpath("./ustb"))
addpath(genpath("./Field_II"))

%% 基本常量
c0 = 1500;      % 声速 [m/s]
fs = 10e6;      % 采样频率 [Hz]
dt = 1/fs;      % 采样步长 [s]

%% Field II 初始化
field_init(0);
set_field('c', c0);              % 声速 [m/s]
set_field('fs', fs);             % 采样频率 [Hz]
set_field('use_rectangles', 1);  % 使用矩形阵元

%% 换能器定义 - Mills Cross配置
%
% MBES使用正交双阵列配置:
% - 发射阵列: 沿跨航向(y轴)放置,形成跨航向宽扇形
% - 接收阵列: 沿沿航向(x轴)放置,用于方位波束形成
%
% 坐标系: x(沿航向-船前进方向), y(跨航向-船左右), z(深度-向下)

% 中心频率和波长
f0 = 200e3;           % 换能器中心频率 [Hz]
lambda = c0/f0;       % 波长 [m]
kerf = 0.05*lambda;   % 阵元间隙 [m]

% 发射阵列参数(跨航向 - 沿y轴排列)
Tx_N = 32;                          % 发射阵元数量
Tx_pitch = lambda;                  % 阵元间距 [m]
Tx_element_width = Tx_pitch - kerf; % 阵元宽度 [m]
Tx_element_height = 15e-3;          % 阵元高度 [m]

% 接收阵列参数(沿航向 - 沿x轴排列)
Rx_N = 128;                          % 接收阵元数量
Rx_pitch = lambda/2;                 % 接收阵间距(更小以提高分辨率) [m]
Rx_element_width = Rx_pitch - kerf;  % 阵元宽度 [m]
Rx_element_height = 15e-3;           % 阵元高度 [m]

disp('========================================');
disp('MBES换能器配置 (Mills Cross)');
disp('========================================');
fprintf('中心频率: %.1f kHz\n', f0/1e3);
fprintf('波长: %.2f mm\n', lambda*1e3);
fprintf('发射阵列(跨航向y轴): %d 阵元, 长度 %.2f cm\n', Tx_N, Tx_N*Tx_pitch*1e2);
fprintf('接收阵列(沿航向x轴): %d 阵元, 长度 %.2f cm\n', Rx_N, Rx_N*Rx_pitch*1e2);
disp('========================================');

%% 海底地形参数定义
seabed_length_along = 60;      % 沿航向长度 [m]
seabed_width_across = 100;     % 跨航向宽度 [m]
seabed_depth_mean = 50;        % 平均水深 [m]
seabed_roughness = 3;          % 海底起伏幅度 [m]
N_seabed_scatter = 3e4;        % 海底散射点数量

% 跨航向角度范围
across_track_angles = linspace(-60, 60, 31);  % 31个角度
Na = length(across_track_angles);

% 帧数
F = 1;  % 单帧仿真

disp(' ');
disp('========================================');
disp('海底场景参数');
disp('========================================');
fprintf('沿航向范围: %.0f m\n', seabed_length_along);
fprintf('跨航向范围: %.0f m\n', seabed_width_across);
fprintf('平均水深: %.0f m\n', seabed_depth_mean);
fprintf('海底起伏: ±%.1f m\n', seabed_roughness);
fprintf('散射点数量: %.0f\n', N_seabed_scatter);
disp('========================================');

%% 脉冲定义 - 使用chirp(线性调频)信号
fractional_bandwidth = 0.5;
BW = fractional_bandwidth * f0;    % 带宽 [Hz]
pulse_len = 1e-3;                  % 脉冲长度 [s]

% 生成chirp激励信号
t_p = 0:dt:pulse_len;
excitation = chirp(t_p, f0-BW/2, pulse_len, f0+BW/2) .* hamming(numel(t_p))';
excitation = excitation - mean(excitation);    % 去除直流分量

% 冲激响应(理想情况)
impulse_response = 1;
one_way_ir = conv(impulse_response, excitation);

% 匹配滤波器(用于脉冲压缩)
MF = conj(flipud(excitation(:)));
two_way_ir = conv(one_way_ir, MF);
lag = length(two_way_ir)/2 + 1;

% 显示脉冲
figure('Name', 'MBES信号分析');
plot((0:(length(two_way_ir)-1))*dt - lag*dt, two_way_ir); 
hold on; grid on; axis tight
plot((0:(length(two_way_ir)-1))*dt - lag*dt, abs(hilbert(two_way_ir)), 'r', 'LineWidth', 1.5)
plot([0 0], [min(two_way_ir) max(two_way_ir)], 'g--', 'LineWidth', 1.5);
legend('双程脉冲', '包络', '估计延迟');
title('Field II 双程冲激响应');
xlabel('时间 [s]');
ylabel('幅度');

%% 创建换能器阵列
%
% 关键修正: 使用xdc_rectangles创建90度旋转的发射阵列
% 接收阵列使用标准xdc_linear_array (沿x轴)
% 发射阵列手动构建 (沿y轴)

% === 接收阵列 (沿x轴,使用标准函数) ===
noSubAz_Rx = max(1, round(Rx_element_width/(lambda/8)));
noSubEl_Rx = max(1, round(Rx_element_height/(lambda/8)));
Rh = xdc_linear_array(Rx_N, Rx_element_width, Rx_element_height, ...
                      kerf, noSubAz_Rx, noSubEl_Rx, [0 0 Inf]); 

% === 发射阵列 (沿y轴,使用xdc_rectangles手动构建) ===
% xdc_rectangles需要19列的rect矩阵:
% [no, x1,y1,z1, x2,y2,z2, x3,y3,z3, x4,y4,z4, apo, width, height, cx,cy,cz]

rect = zeros(Tx_N, 19);
center_tx = zeros(Tx_N, 3);

% 计算阵元中心位置(沿y轴排列)
y_positions = ((1:Tx_N) - (Tx_N+1)/2) * Tx_pitch;

for i = 1:Tx_N
    % 阵元编号
    rect(i, 1) = i;
    
    % 阵元中心
    cx = 0;
    cy = y_positions(i);
    cz = 0;
    
    % 阵元半宽半高
    hw = Tx_element_height/2;  % 沿x方向(旋转后)
    hh = Tx_element_width/2;   % 沿y方向(旋转后)
    
    % 四个角点坐标(逆时针顺序,从俯视图看)
    % 角点1: (-hw, cy-hh, 0)
    % 角点2: (+hw, cy-hh, 0)
    % 角点3: (+hw, cy+hh, 0)
    % 角点4: (-hw, cy+hh, 0)
    rect(i, 2:4)   = [-hw, cy-hh, 0];  % 角点1
    rect(i, 5:7)   = [+hw, cy-hh, 0];  % 角点2
    rect(i, 8:10)  = [+hw, cy+hh, 0];  % 角点3
    rect(i, 11:13) = [-hw, cy+hh, 0];  % 角点4
    
    % 加权值
    rect(i, 14) = 1.0;
    
    % 阵元尺寸
    rect(i, 15) = Tx_element_height;  % width (x方向)
    rect(i, 16) = Tx_element_width;   % height (y方向)
    
    % 中心点
    rect(i, 17:19) = [cx, cy, cz];
    
    % 物理阵元中心(用于center参数)
    center_tx(i, :) = [cx, cy, cz];
end

% 创建发射阵列
Th = xdc_rectangles(rect, center_tx, [0 0 Inf]);

% 设置激励和冲激响应
xdc_excitation(Th, excitation);
xdc_impulse(Th, impulse_response);
xdc_baffle(Th, 0);
xdc_center_focus(Th, [0 0 0]);

xdc_impulse(Rh, impulse_response);
xdc_baffle(Rh, 0);
xdc_center_focus(Rh, [0 0 0]);

disp('换能器阵列创建完成');

%% 海底地形生成
disp(' ');
disp('正在生成海底地形...');

% 计算网格点数
Nx_seabed = round(sqrt(N_seabed_scatter * seabed_length_along / seabed_width_across));
Ny_seabed = round(N_seabed_scatter / Nx_seabed);

% 生成均匀网格
x_seabed = linspace(-seabed_length_along/2, seabed_length_along/2, Nx_seabed);
y_seabed = linspace(-seabed_width_across/2, seabed_width_across/2, Ny_seabed);
[Xg, Yg] = meshgrid(x_seabed, y_seabed);

% 生成真实的海底地形(多尺度起伏)
Zg = seabed_depth_mean + ...
     seabed_roughness * (1.0*sin(2*pi*Xg/30) + ...
                         1.2*sin(2*pi*Yg/40) + ...
                         0.8*sin(2*pi*Xg/10) .* cos(2*pi*Yg/15) + ...
                         0.5*randn(size(Xg)));

% 添加一个小山丘特征
hill_x = 0;
hill_y = 0;
hill_height = 15;
hill_radius = 15;
distance = sqrt((Xg-hill_x).^2 + (Yg-hill_y).^2);
Zg = Zg - hill_height * exp(-(distance/hill_radius).^2);

% 展开为点列表 - 关键修正: 定义point_position和point_amplitudes
point_position = [Xg(:), Yg(:), Zg(:)];  % [沿航向, 跨航向, 深度]

% 设置海底散射系数
point_amplitudes = 0.1 + 0.2 * abs(randn(size(point_position,1), 1));
point_amplitudes = point_amplitudes .* (1 + 0.3*sin(2*pi*Yg(:)/20));

disp(['海底地形生成完成: ', num2str(size(point_position,1)), ' 个散射点']);

%% 场景可视化
disp('正在生成场景预览...');

figure('Name', 'MBES场景预览', 'Position', [50 60 1400 700]);

% 子图1: 3D视图
subplot(2,2,1);
scatter3(point_position(1:20:end,1), point_position(1:20:end,2), point_position(1:20:end,3), ...
         3, point_position(1:20:end,3), 'filled'); 
hold on;
plot3(0, 0, 0, 'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'r');
% 发射阵列方向(y轴)
quiver3(0, 0, 0, 0, 10, 0, 'r', 'LineWidth', 2, 'MaxHeadSize', 1);
% 接收阵列方向(x轴)
quiver3(0, 0, 0, 10, 0, 0, 'b', 'LineWidth', 2, 'MaxHeadSize', 1);
set(gca, 'ZDir', 'reverse');
axis equal tight; grid on;
xlabel('沿航向 x [m]'); ylabel('跨航向 y [m]'); zlabel('深度 z [m]');
title('3D场景视图');
view(35, 25);
colorbar;
legend('海底', '换能器', 'Tx(y轴)', 'Rx(x轴)', 'Location', 'best');

% 子图2: 顶视图
subplot(2,2,2);
scatter(point_position(1:20:end,1), point_position(1:20:end,2), ...
        3, point_position(1:20:end,3), 'filled');
hold on;
plot(0, 0, 'rp', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
% 绘制波束覆盖范围
swath_angle = 60;
y_coverage = seabed_depth_mean * tand(swath_angle);
plot([0, 0], [-y_coverage, y_coverage], 'r--', 'LineWidth', 2);
axis equal tight; grid on;
xlabel('沿航向 x [m]'); ylabel('跨航向 y [m]');
title('顶视图 (x-y)');
colorbar;

% 子图3: 跨航向剖面(x≈0)
subplot(2,2,3);
idx_profile = abs(Xg(:)) < 2;
scatter(point_position(idx_profile,2), point_position(idx_profile,3), ...
        10, point_position(idx_profile,3), 'filled');
hold on;
plot(0, 0, 'rp', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
set(gca, 'YDir', 'reverse');
axis equal tight; grid on;
xlabel('跨航向 y [m]'); ylabel('深度 z [m]');
title('跨航向剖面 (x≈0)');
colorbar;

% 子图4: 换能器阵列示意图
subplot(2,2,4);
hold on; grid on;
% 绘制发射阵元(沿y轴)
for i = 1:Tx_N
    rectangle('Position', [-Tx_element_height/2*100, (y_positions(i)-Tx_element_width/2)*100, ...
              Tx_element_height*100, Tx_element_width*100], ...
              'FaceColor', [1 0.5 0.5], 'EdgeColor', 'r');
end
% 绘制接收阵元(沿x轴)
x_rx_positions = ((1:Rx_N) - (Rx_N+1)/2) * Rx_pitch;
for i = 1:Rx_N
    rectangle('Position', [(x_rx_positions(i)-Rx_element_width/2)*100, -Rx_element_height/2*100, ...
              Rx_element_width*100, Rx_element_height*100], ...
              'FaceColor', [0.5 0.5 1], 'EdgeColor', 'b');
end
axis equal;
xlabel('x [cm]'); ylabel('y [cm]');
title('Mills Cross阵列布局');
legend({'Tx阵列(y轴)', 'Rx阵列(x轴)'}, 'Location', 'best');
xlim([-30 30]); ylim([-30 30]);

drawnow;

%% Field II 仿真计算
disp(' ');
disp('========================================');
disp('开始Field II仿真计算...');
disp('========================================');

time_start = tic;
wb = waitbar(0, 'MBES: 正在计算跨航向波束数据...');

% Hamming窗用于接收孔径加权
win = hamming(Rx_N);

% 深度轴定义
depth_axis = linspace(seabed_depth_mean-20, seabed_depth_mean+20, 256)';
n_depth = length(depth_axis);

% 初始化波束形成图像
bf_image = zeros(n_depth, Na);

% 记录海底检测数据
max_amplitudes = zeros(Na, 1);
bottom_times = zeros(Na, 1);

for f = 1:F
    waitbar(0, wb, sprintf('MBES: 正在计算数据集, 第 %d/%d 帧', f, F));
    
    for n = 1:Na
        waitbar(n/Na, wb);
        
        angle = across_track_angles(n);
        
        % === 发射阵列设置 ===
        % 所有阵元同时发射,不聚焦(平面波)
        xdc_apodization(Th, 0, ones(1, Tx_N));
        
        % 发射聚焦到无穷远(平面波)或指定方向
        % 对于跨航向角度angle,聚焦点在y-z平面内
        focus_depth = 1000;  % 远场聚焦
        focus_x = 0;
        focus_y = focus_depth * sind(angle);
        focus_z = focus_depth * cosd(angle);
        xdc_focus(Th, 0, [focus_x, focus_y, focus_z]);
        
        % === 接收阵列设置 ===
        xdc_apodization(Rh, 0, ones(1, Rx_N));
        xdc_focus_times(Rh, 0, zeros(1, Rx_N));  % 动态聚焦
        
        % === 执行散射计算 ===
        [v, t] = calc_scat_multi(Th, Rh, point_position, point_amplitudes);
        
        % 检查数据有效性
        if isempty(v) || all(v(:) == 0)
            fprintf('  警告: 角度 %.1f° 无有效回波\n', angle);
            continue;
        end
        
        % === 匹配滤波 + 加窗 ===
        sig_mf = zeros(size(v));
        for m = 1:Rx_N
            tmp_sig = conv(v(:,m), MF, 'same');
            sig_mf(:,m) = tmp_sig * win(m);
        end
        
        % === 波束形成(DAS) ===
        beam_signal = sum(sig_mf, 2);
        
        % 时间轴
        time_axis = (0:size(beam_signal,1)-1) * dt + t - lag*dt;
        
        % 转换为斜距
        range_axis = time_axis * c0 / 2;
        
        % 斜距转换为垂直深度(考虑角度)
        depth_from_range = range_axis * cosd(angle);
        
        % 插值到统一深度轴
        beam_signal_interp = interp1(depth_from_range, abs(beam_signal), ...
                                     depth_axis, 'linear', 0);
        
        % 存储波束形成结果
        bf_image(:, n) = beam_signal_interp;
        
        % 记录最大回波位置
        [max_amp, max_idx] = max(abs(beam_signal));
        max_amplitudes(n) = max_amp;
        if max_idx <= length(range_axis) && max_idx > 0
            bottom_times(n) = range_axis(max_idx);
        else
            bottom_times(n) = NaN;
        end
        
        % 清理临时变量
        clear v sig_mf tmp_sig beam_signal;
        
        % 进度显示
        if mod(n, 5) == 0 || n == Na
            fprintf('  完成: %d/%d (%.1f%%)\n', n, Na, n/Na*100);
        end
    end
end
close(wb);

time_elapsed = toc(time_start);
disp(['Field II计算完成! 用时: ', num2str(time_elapsed, '%.1f'), ' 秒']);

% 归一化波束形成图像
if max(bf_image(:)) > 0
    bf_image = bf_image / max(bf_image(:));
end

disp(' ');
disp('========================================');
disp('波束形成数据统计');
disp('========================================');
fprintf('数据维度: [深度, 角度] = [%d, %d]\n', n_depth, Na);
fprintf('内存占用: %.2f MB\n', numel(bf_image)*8/1e6);
disp('========================================');

%% 成像显示
% 转换为dB
bf_image_db = 20*log10(abs(bf_image) + eps);
bf_image_db = bf_image_db - max(bf_image_db(:));
bf_image_db(bf_image_db < -60) = -60;

disp(' ');
disp('正在生成成像结果...');

figure('Name', 'MBES成像结果', 'Position', [100, 50, 1400, 700]);

% 子图1: 角度-深度图
subplot(1,2,1);
imagesc(across_track_angles, depth_axis, bf_image_db);
colormap(hot);
colorbar;
clim([-60 0]);
xlabel('跨航向角度 [°]');
ylabel('深度 [m]');
title('MBES波束形成图像 (角度-深度)');
set(gca, 'YDir', 'reverse');
grid on;

% 子图2: 直角坐标
subplot(1,2,2);
[ANGLE_grid, DEPTH_grid] = meshgrid(across_track_angles*pi/180, depth_axis);
Y_grid = DEPTH_grid .* sin(ANGLE_grid);
Z_grid = DEPTH_grid .* cos(ANGLE_grid);

pcolor(Y_grid, Z_grid, bf_image_db);
shading interp;
colormap(hot);
colorbar;
clim([-60 0]);
xlabel('跨航向距离 [m]');
ylabel('深度 [m]');
title('MBES成像结果 (直角坐标)');
axis equal tight;
set(gca, 'YDir', 'reverse');
grid on;

%% 3D海底重建
disp('正在进行3D海底重建...');

figure('Name', 'MBES 3D海底重建', 'Position', [100, 100, 1400, 700]);

% 子图1: 波束-深度图
subplot(2,2,1);
imagesc(across_track_angles, depth_axis, bf_image_db);
colormap(hot);
colorbar;
clim([-60 0]);
xlabel('跨航向角度 [°]');
ylabel('深度 [m]');
title('波束-深度图');
set(gca, 'YDir', 'reverse');
grid on;

% 子图2: 海底检测
subplot(2,2,2);
bottom_depths = bottom_times;
valid_idx = ~isnan(bottom_depths) & bottom_depths > 0;

plot(across_track_angles(valid_idx), bottom_depths(valid_idx), ...
     'b.-', 'LineWidth', 1.5, 'MarkerSize', 8);
hold on;
plot(across_track_angles, seabed_depth_mean*ones(size(across_track_angles)), ...
     'r--', 'LineWidth', 2);
grid on;
xlabel('跨航向角度 [°]');
ylabel('检测深度 [m]');
title('海底深度检测');
legend('检测深度', '平均深度', 'Location', 'best');
set(gca, 'YDir', 'reverse');
axis tight;

% 子图3: 3D点云
subplot(2,2,3);
angles_valid = across_track_angles(valid_idx);
depths_valid = bottom_depths(valid_idx);

x_bottom = zeros(size(depths_valid(:)));
y_bottom = depths_valid(:) .* sind(angles_valid(:));
z_bottom = depths_valid(:) .* cosd(angles_valid(:));

scatter3(x_bottom, y_bottom, z_bottom, 30, 'b', 'filled');
hold on;
plot3(0, 0, 0, 'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'r');

xlabel('沿航向 x [m]');
ylabel('跨航向 y [m]');
zlabel('深度 z [m]');
title('3D测深点云');
set(gca, 'ZDir', 'reverse');
axis equal;
grid on;
view(35, 25);

% 子图4: 与真实海底对比
subplot(2,2,4);
idx_true = abs(Xg(:)) < 2;
plot(point_position(idx_true, 2), point_position(idx_true, 3), 'k.', 'MarkerSize', 3);
hold on;
plot(y_bottom, z_bottom, 'b.-', 'LineWidth', 2, 'MarkerSize', 10);
xlabel('跨航向 y [m]');
ylabel('深度 z [m]');
title('海底重建对比');
legend('真实海底', '重建海底', 'Location', 'best');
set(gca, 'YDir', 'reverse');
grid on;
axis tight;

%% 统计分析
disp(' ');
disp('========================================');
disp('成像质量评估');
disp('========================================');

valid_idx = ~isnan(bottom_depths) & bottom_depths > 0;
angles_valid = across_track_angles(valid_idx) * pi/180;
depths_valid = bottom_depths(valid_idx);

y_bottom_stat = depths_valid .* sin(angles_valid);
swath_width = max(y_bottom_stat) - min(y_bottom_stat);
coverage_ratio = swath_width / seabed_depth_mean;

fprintf('测线覆盖宽度: %.1f m\n', swath_width);
fprintf('覆盖比(宽度/水深): %.2f\n', coverage_ratio);
fprintf('有效测深点数: %d / %d (%.1f%%)\n', ...
        sum(valid_idx), length(valid_idx), sum(valid_idx)/length(valid_idx)*100);
fprintf('总用时: %.1f 秒\n', time_elapsed);

disp('========================================');

%% 释放Field II资源
field_end;

disp(' ');
disp('========================================');
disp('MBES仿真处理完成!');
disp('========================================');