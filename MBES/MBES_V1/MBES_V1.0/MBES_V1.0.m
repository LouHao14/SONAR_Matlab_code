%% 使用Field II和USTB计算多波束测深声纳(MBES)数据集并进行波束形成
%
% 本示例展示如何使用Field II创建多波束测深声纳(MBES)成像仿真,
% 并将数据存储到USTB的channel_data对象中,然后使用波束形成算法。
% 
% 主要功能:
% - Mills Cross正交阵列配置(发射阵+接收阵)
% - 生成海底地形散射体场景
% - 使用chirp信号进行声场仿真
% - 俯仰向扇形扫描
% - 多波束形成算法
% - 3D点云和地形网格显示
%
% MBES特点:
% - 发射阵列: 跨航向(y轴),形成俯仰向窄波束
% - 接收阵列: 沿航向(x轴),形成多个方位波束
% - 输出: 离散的角度-距离测量点,用于海底地形测绘
%
% Field II仿真程序(field-ii.dk)需要在MATLAB路径中
%
% 日期:     2025
% 作者:     [您的名字]

%% 清除工作空间和关闭图形窗口

clear all;
close all;

%% 添加必要的工具箱路径

addpath(genpath("./ustb"))
addpath(genpath("./Field_II"))

%% 基本常量
% 
% 首先定义成像场景的基本常量,包括声速、采样频率和采样步长

c0 = 1500;      % 声速 [m/s]
fs = 10e6;      % 采样频率 [Hz]
dt = 1/fs;      % 采样步长 [s]

%% Field II 初始化
% 
% 初始化Field II工具箱。这需要Field II仿真程序(<field-ii.dk>)在MATLAB路径中。
% 同时将设定的常量传递给它。

field_init(0);
set_field('c',c0);              % 声速 [m/s]
set_field('fs',fs);             % 采样频率 [Hz]
set_field('use_rectangles',1);  % 使用矩形阵元

%% 换能器定义 - Mills Cross 正交双阵列配置
% 
% MBES使用两个正交的线阵:
% - 发射阵列: 沿跨航向(y轴),用于俯仰向聚焦
% - 接收阵列: 沿航向(x轴),用于方位向波束形成
% 这种配置可以形成窄的"扇形"波束覆盖

% 基本参数
f0 = 0.9e+06;          % 换能器中心频率 [Hz]
lambda = c0/f0;        % 波长 [m]
kerf = 0.1*lambda;     % 阵元间隙 [m]

% 发射阵列 (跨航向, y轴)
Tx_probe = uff.linear_array();
Tx_probe.N = 32;                           % 发射阵元数
Tx_probe.pitch = lambda;                   % 阵元间距 [m]
Tx_probe.element_width = Tx_probe.pitch - kerf;   % 阵元宽度 [m]
Tx_probe.element_height = 15e-3;           % 阵元高度 [m]

% 接收阵列 (沿航向, x轴)
Rx_probe = uff.linear_array();
Rx_probe.N = 128;                          % 接收阵元数(更多,用于高角度分辨率)
Rx_probe.pitch = lambda/2;                 % 接收阵间距更小 [m]
Rx_probe.element_width = Rx_probe.pitch - kerf;   % 阵元宽度 [m]
Rx_probe.element_height = 15e-3;           % 阵元高度 [m]

pulse_duration = 2.5;                      % 脉冲持续时间 [周期]

disp('=== 多波束测深声纳配置 ===');
disp(['发射阵列: ', num2str(Tx_probe.N), ' 阵元 (跨航向)']);
disp(['接收阵列: ', num2str(Rx_probe.N), ' 阵元 (沿航向)']);
disp(['中心频率: ', num2str(f0/1e6), ' MHz']);

%% 测量参数定义

% 俯仰角扫描范围(发射)
pitch_angles = -60:2:60;           % 俯仰角: -60° 到 +60°, 步进2° [度]
N_pitch = length(pitch_angles);    % 俯仰发射数量

% 方位角范围(接收波束形成)
azimuth_angles = -60:1:60;         % 方位角: -60° 到 +60°, 步进1° [度]
N_azimuth = length(azimuth_angles);% 方位波束数量

% 测深范围
max_range = 100;                   % 最大测量距离 [m]
focus_range = 50;                  % 发射聚焦距离 [m]

disp(['俯仰扫描: ', num2str(N_pitch), ' 个角度']);
disp(['方位波束: ', num2str(N_azimuth), ' 个波束']);
disp(['测深范围: 0-', num2str(max_range), ' m']);

%% 海底地形参数定义
% 
% 生成真实的海底地形用于测深仿真

seabed_len = 100;          % 海底长度(沿航向, x方向) [m]
seabed_wid = 80;           % 海底宽度(跨航向, y方向) [m]
seabed_depth0 = 50;        % 平均水深 [m]
seabed_amp = 5;            % 海底起伏幅度 [m]
N_seabed_scatter = 5e4;    % 海底散射点数量

disp(['海底尺寸: ', num2str(seabed_len), 'm × ', num2str(seabed_wid), 'm']);
disp(['平均水深: ', num2str(seabed_depth0), ' m']);

%% 脉冲定义
% 
% 使用chirp(线性调频)信号定义脉冲-回波信号。

pulse = uff.pulse();
pulse.fractional_bandwidth = 1/9;          % 相对带宽
pulse.center_frequency = f0;               % 中心频率 [Hz]
BW = pulse.fractional_bandwidth * f0;      % 带宽 [Hz]
pulse_len = 2e-3;                          % 脉冲长度 [s]

% 生成chirp激励信号
t_p = 0:dt:pulse_len;
excitation = chirp(t_p, f0-BW/2, pulse_len, f0+BW/2) .* hamming(numel(t_p))';
excitation = excitation - mean(excitation);  % 去除直流分量

% 冲激响应(理想情况)
impulse_response = 1;
one_way_ir = conv(impulse_response, excitation);

% 匹配滤波器(用于脉冲压缩)
MF = conj(flipud(excitation(:)));
two_way_ir = conv(one_way_ir, MF);
lag = length(two_way_ir)/2 + 1;

% 显示脉冲
figure('Name', 'MBES脉冲特性');
plot((0:(length(two_way_ir)-1))*dt - lag*dt, two_way_ir); 
hold on; grid on; axis tight
plot((0:(length(two_way_ir)-1))*dt - lag*dt, abs(hilbert(two_way_ir)), 'r', 'LineWidth', 1.5)
plot([0 0], [min(two_way_ir) max(two_way_ir)], 'g--', 'LineWidth', 2);
legend('双程脉冲', '包络', '估计延迟');
xlabel('时间 [s]'); ylabel('幅度');
title('Field II 双程冲激响应');

%% 孔径对象 - 创建发射和接收阵列
% 
% 使用Field II创建两个阵列
% 简化方案: 使用两个标准线阵,通过聚焦策略实现MBES效果
% 
% 注意: 由于Field II某些版本不支持阵列旋转函数,
%       这里采用简化的配置方式

% 计算子阵元数
noSubAz_Tx = round(Tx_probe.element_width/(lambda/8));
noSubEl_Tx = round(Tx_probe.element_height/(lambda/8));
noSubAz_Rx = round(Rx_probe.element_width/(lambda/8));
noSubEl_Rx = round(Rx_probe.element_height/(lambda/8));

% === 发射阵列 (用于俯仰向聚焦) ===
% 创建较短的发射阵列
Th = xdc_linear_array(Tx_probe.N, Tx_probe.element_width, Tx_probe.element_height, ...
                      kerf, noSubAz_Tx, noSubEl_Tx, [0 0 Inf]);

% === 接收阵列 (用于方位向波束形成) ===
% 创建较长的接收阵列
Rh = xdc_linear_array(Rx_probe.N, Rx_probe.element_width, Rx_probe.element_height, ...
                      kerf, noSubAz_Rx, noSubEl_Rx, [0 0 Inf]); 

% 设置激励、冲激响应和障板
xdc_excitation(Th, excitation);
xdc_impulse(Th, impulse_response);
xdc_baffle(Th, 0);
xdc_center_focus(Th, [0 0 0]);

xdc_impulse(Rh, impulse_response);
xdc_baffle(Rh, 0);
xdc_center_focus(Rh, [0 0 0]);

disp('阵列创建完成 - 简化MBES配置');
disp(['  发射阵列: ', num2str(Tx_probe.N), ' 阵元']);
disp(['  接收阵列: ', num2str(Rx_probe.N), ' 阵元']);
disp('  注意: 使用简化配置,两阵列平行放置');
disp('        通过发射聚焦和接收波束形成实现测深效果');

%% 定义海底散射体场景
% 
% 生成具有真实起伏的海底地形散射点

disp('正在生成海底地形...');

% 计算网格点数
Nx = round(sqrt(N_seabed_scatter * seabed_len / seabed_wid));
Ny = round(N_seabed_scatter / Nx);

% 生成均匀网格 (x: 沿航向, y: 跨航向)
x_lin = linspace(0, seabed_len, Nx);           % 沿航向 0-100m
y_lin = linspace(-seabed_wid/2, seabed_wid/2, Ny); % 跨航向 -40 到 +40m
[Xg, Yg] = meshgrid(x_lin, y_lin);

% 生成海底深度起伏
% 多种地形特征的组合
Zg = seabed_depth0 + ...
     seabed_amp * 0.3 * sin(2*pi*Xg/seabed_len) + ...      % 沿航向波动
     seabed_amp * 0.3 * sin(2*pi*Yg/seabed_wid) + ...      % 跨航向波动
     seabed_amp * 0.2 * sin(4*pi*Xg/seabed_len) .* cos(2*pi*Yg/seabed_wid) + ... % 交叉波纹
     seabed_amp * 0.2 * randn(size(Xg));                    % 随机粗糙度

% 添加一个海底山丘特征
hill_x = seabed_len/2;
hill_y = 0;
hill_amp = 8;
hill_width = 15;
dist_to_hill = sqrt((Xg - hill_x).^2 + (Yg - hill_y).^2);
Zg = Zg - hill_amp * exp(-(dist_to_hill/hill_width).^2);

% 展开为点列表 [x, y, z] (z为深度,向下为正)
pos_seabed = [Xg(:) Yg(:) Zg(:)];

% 设置海底散射系数 (随机变化模拟不同底质)
amp_seabed = 0.15 + 0.1 * abs(randn(size(pos_seabed,1), 1));

disp(['海底散射体生成完成: ', num2str(size(pos_seabed,1)), ' 个散射点']);

%% 场景可视化
% 
% 显示完整的3D海底地形场景

disp('正在生成场景预览...');

figure('Name', '多波束测深声纳场景预览', 'Position', [50 60 1400 700]);

% 子图1: 3D海底地形
subplot(2,2,1);
% 使用surf显示地形(降采样以提高显示速度)
downsample_factor = 5;
surf(Xg(1:downsample_factor:end, 1:downsample_factor:end), ...
     Yg(1:downsample_factor:end, 1:downsample_factor:end), ...
     Zg(1:downsample_factor:end, 1:downsample_factor:end));
hold on;
plot3(0, 0, 0, 'rp', 'MarkerFaceColor', 'r', 'MarkerSize', 15);
shading interp;
colormap(jet);
colorbar;
set(gca, 'ZDir', 'reverse');
axis tight; grid on;
xlabel('沿航向 x [m]'); ylabel('跨航向 y [m]'); zlabel('深度 z [m]');
title('3D 海底地形');
view(45, 30);

% 子图2: 顶视图
subplot(2,2,2);
contourf(Xg, Yg, Zg, 20);
hold on;
plot(0, 0, 'rp', 'MarkerFaceColor', 'r', 'MarkerSize', 12);
colorbar;
axis equal tight; grid on;
xlabel('沿航向 x [m]'); ylabel('跨航向 y [m]');
title('顶视图 - 深度等值线');

% 子图3: 沿航向剖面
subplot(2,2,3);
plot(x_lin, Zg(round(Ny/2), :), 'b-', 'LineWidth', 1.5);
hold on;
plot([0 seabed_len], [seabed_depth0 seabed_depth0], 'r--', 'LineWidth', 1);
set(gca, 'YDir', 'reverse');
grid on;
xlabel('沿航向 x [m]'); ylabel('深度 z [m]');
title('沿航向中线剖面');
legend('海底地形', '平均深度');

% 子图4: 跨航向剖面
subplot(2,2,4);
plot(y_lin, Zg(:, round(Nx/2)), 'b-', 'LineWidth', 1.5);
hold on;
plot([-seabed_wid/2 seabed_wid/2], [seabed_depth0 seabed_depth0], 'r--', 'LineWidth', 1);
set(gca, 'YDir', 'reverse');
grid on;
xlabel('跨航向 y [m]'); ylabel('深度 z [m]');
title('跨航向中线剖面');
legend('海底地形', '平均深度');

drawnow;

%% 显示波束覆盖模式
% 
% 可视化MBES的俯仰-方位波束覆盖

figure('Name', 'MBES波束覆盖模式', 'Position', [100 100 1200 500]);

% 子图1: 俯仰角覆盖(侧视图)
subplot(1,2,1);
hold on; grid on;
for p = 1:5:N_pitch  % 每隔5个角度绘制一条
    angle = pitch_angles(p) * pi/180;
    x_ray = [0, focus_range * cos(angle)];
    z_ray = [0, focus_range * sin(angle)];
    plot(x_ray, z_ray, 'b-', 'LineWidth', 0.5);
end
plot(0, 0, 'rp', 'MarkerFaceColor', 'r', 'MarkerSize', 12);
set(gca, 'YDir', 'reverse');
axis equal;
xlim([0 max_range]); ylim([-max_range/2 max_range/2]);
xlabel('距离 [m]'); ylabel('深度变化 [m]');
title('俯仰角覆盖 (沿航向视图)');

% 子图2: 方位角覆盖(顶视图)
subplot(1,2,2);
hold on; grid on;
for a = 1:10:N_azimuth  % 每隔10个角度绘制一条
    angle = azimuth_angles(a) * pi/180;
    x_ray = [0, max_range * cos(angle)];
    y_ray = [0, max_range * sin(angle)];
    plot(x_ray, y_ray, 'g-', 'LineWidth', 0.5);
end
plot(0, 0, 'rp', 'MarkerFaceColor', 'r', 'MarkerSize', 12);
axis equal;
xlim([0 max_range]); ylim([-max_range max_range]);
xlabel('沿航向 x [m]'); ylabel('跨航向 y [m]');
title('方位角覆盖 (顶视图)');

drawnow;

%% 组合所有散射点

point_position = pos_seabed;
point_amplitudes = amp_seabed;

%% 输出数据初始化

% 计算合理的时间窗长度
% 注意: 避免过大的cropat导致内存溢出
max_time = 2 * max_range / c0;  % 双程传播时间
cropat_ideal = round(max_time / dt);   % 理想采样点数

% 根据测深范围设置合理的限制
% cropat的对应关系: 1000采样点 ≈ 7.5m双程距离
if max_range <= 50
    max_samples = 15000;  % 50m对应约13,333采样点
elseif max_range <= 100
    max_samples = 30000;  % 100m对应约26,667采样点
else
    max_samples = 50000;  % 更大范围
end

cropat = min(cropat_ideal, max_samples);

% 估算内存使用
memory_mb = cropat * Rx_probe.N * N_pitch * 8 / 1e6;  % 8 bytes per double
disp(['cropat: ', num2str(cropat), ' 采样点 (对应约 ', ...
      num2str(cropat*dt*c0/2, '%.1f'), ' m)']);
disp(['CPW内存: ', num2str(memory_mb, '%.1f'), ' MB']);

if memory_mb > 2000
    warning('CPW数组超过2GB, 可能导致内存不足');
    disp('建议: 减少N_pitch或max_range参数');
elseif memory_mb > 1000
    warning('CPW数组超过1GB, 请确保有足够内存');
end

CPW = zeros(cropat, Rx_probe.N, N_pitch, 1);  % 通道数据矩阵

%% 计算多波束通道数据
% 
% 对每个俯仰角进行发射,记录所有接收通道的数据

disp(' ');
disp('========================================');
disp('开始Field II多波束数据计算');
disp('========================================');

wb = waitbar(0, 'Field II: 正在计算多波束测深数据集');

% Blackman窗用于接收孔径加权
win = blackman(Rx_probe.N);

tic;  % 开始计时

for p = 1:N_pitch
    waitbar(p/N_pitch, wb, sprintf('Field II: 计算俯仰角 %d/%d (%.1f°)', ...
                                    p, N_pitch, pitch_angles(p)));
    
    % 当前俯仰角
    pitch_rad = pitch_angles(p) * pi/180;
    
    % 计算聚焦点坐标
    % 在简化配置中,我们使用发射阵列的横向聚焦来模拟俯仰效果
    % 聚焦点: 在z轴方向偏移,模拟不同的俯仰角
    focus_x = 0;
    focus_y = 0;  
    focus_z = focus_range;  % 主要聚焦距离
    
    % 使用横向偏移模拟俯仰角效果
    % 通过调整发射阵元的激活模式和延迟
    focus_point = [focus_x, focus_y, focus_z];
    
    % 发射孔径设置
    % 根据俯仰角调整发射阵元的加权
    tx_weights = ones(1, Tx_probe.N);
    % 可以添加一个渐变加权来改善指向性
    center_idx = (Tx_probe.N + 1) / 2;
    for i = 1:Tx_probe.N
        offset = (i - center_idx) * Tx_probe.pitch;
        % 根据俯仰角计算相位延迟
        phase_delay = offset * sin(pitch_rad) / c0;
        tx_weights(i) = 1.0;  % 保持均匀加权,通过延迟控制方向
    end
    
    xdc_apodization(Th, 0, tx_weights);
    
    % 计算发射聚焦延迟 - 使其指向特定俯仰角
    tx_elem_positions = ((0:Tx_probe.N-1) - (Tx_probe.N-1)/2) * Tx_probe.pitch;
    tx_focus_delays = tx_elem_positions * sin(pitch_rad) / c0;
    % 转换为相对延迟(最大延迟为0)
    tx_focus_delays = tx_focus_delays - min(tx_focus_delays);
    xdc_times_focus(Th, 0, tx_focus_delays);
    
    % 接收孔径设置 - 全孔径接收,不聚焦
    xdc_apodization(Rh, 0, ones(1, Rx_probe.N));
    xdc_focus_times(Rh, 0, zeros(1, Rx_probe.N));
    
    % 执行散射计算
    [v, t] = calc_scat_multi(Th, Rh, point_position, point_amplitudes);
    
    % 重要: 检查返回数据大小，Field II可能返回过大的数组
    if size(v, 1) > cropat
        % 如果返回数据超过cropat，只取前cropat个点
        v = v(1:cropat, :);
    end
    
    % 匹配滤波 + 加窗 (降低旁瓣)
    sig_mf = zeros(size(v));
    for m = 1:Rx_probe.N
        tmp_sig = conv(v(:,m), MF, 'same');
        sig_mf(:,m) = tmp_sig * win(m);
    end
    
    % 存储数据 - 确保不超过cropat
    data_len = min(size(sig_mf, 1), cropat);
    CPW(1:data_len, :, p, 1) = sig_mf(1:data_len, :);
    
    % 保存发射序列信息
    seq(p) = uff.wave();
    seq(p).probe = Rx_probe;  % 使用接收阵列作为参考
    seq(p).source.azimuth = 0;           % 方位角(沿航向方向)
    seq(p).source.elevation = pitch_rad;  % 俯仰角
    seq(p).source.distance = Inf;         % 近似平面波
    seq(p).sound_speed = c0;
    seq(p).delay = -lag*dt + t;
end

close(wb);

elapsed_time = toc;
disp(['数据计算完成! 耗时: ', num2str(elapsed_time, '%.1f'), ' 秒']);

%% 通道数据存储
% 
% 创建UFF数据结构来存储捕获的超声通道数据

channel_data = uff.channel_data();
channel_data.sampling_frequency = fs;
channel_data.sound_speed = c0;
channel_data.initial_time = 0;
channel_data.pulse = pulse;
channel_data.probe = Rx_probe;
channel_data.sequence = seq;
channel_data.data = CPW ./ max(CPW(:));  % 归一化

disp('通道数据已存储到UFF格式');

%% 多波束波束形成
% 
% 对每个俯仰发射,形成多个方位接收波束
% 使用延迟求和(DAS)算法

disp(' ');
disp('========================================');
disp('开始多波束波束形成');
disp('========================================');

% 初始化波束数据
beam_data = zeros(N_pitch, N_azimuth, cropat);

wb = waitbar(0, '正在进行波束形成...');

for p = 1:N_pitch
    waitbar(p/N_pitch, wb, sprintf('波束形成: 俯仰角 %d/%d', p, N_pitch));
    
    for a = 1:N_azimuth
        azimuth_rad = azimuth_angles(a) * pi/180;
        
        % 计算方位向的相位延迟
        % 接收阵列沿x轴,阵元位置
        element_positions = (0:Rx_probe.N-1) * Rx_probe.pitch - ...
                           (Rx_probe.N-1) * Rx_probe.pitch / 2;
        
        % 计算到达时间差
        time_delays = element_positions * sin(azimuth_rad) / c0;
        sample_delays = round(time_delays / dt);
        
        % 延迟求和波束形成
        beam_signal = zeros(cropat, 1);
        for m = 1:Rx_probe.N
            delay = sample_delays(m);
            if delay >= 0
                shifted_signal = [zeros(delay, 1); CPW(1:cropat-delay, m, p, 1)];
            else
                shifted_signal = [CPW(-delay+1:cropat, m, p, 1); zeros(-delay, 1)];
            end
            beam_signal = beam_signal + shifted_signal * win(m);
        end
        
        % 存储波束数据
        beam_data(p, a, :) = beam_signal;
    end
end

close(wb);
disp('波束形成完成!');

%% 海底检测
% 
% 从每个波束中检测海底回波,提取角度-距离-深度信息

disp(' ');
disp('========================================');
disp('开始海底回波检测');
disp('========================================');

% 初始化检测结果
sounding_data = struct();
sounding_data.pitch = [];
sounding_data.azimuth = [];
sounding_data.range = [];
sounding_data.amplitude = [];
sounding_data.x = [];
sounding_data.y = [];
sounding_data.z = [];

% 时间轴
time_axis = (0:cropat-1) * dt;

wb = waitbar(0, '正在检测海底回波...');

detection_count = 0;

for p = 1:N_pitch
    waitbar(p/N_pitch, wb);
    
    pitch_rad = pitch_angles(p) * pi/180;
    
    for a = 1:N_azimuth
        azimuth_rad = azimuth_angles(a) * pi/180;
        
        % 获取波束信号包络
        beam_signal = squeeze(beam_data(p, a, :));
        envelope = abs(hilbert(beam_signal));
        
        % 简单阈值检测 (寻找第一个强回波)
        threshold = 0.3 * max(envelope);
        detection_indices = find(envelope > threshold);
        
        if ~isempty(detection_indices)
            % 取第一个检测点作为海底
            detect_idx = detection_indices(1);
            detect_time = time_axis(detect_idx);
            detect_range = detect_time * c0 / 2;  % 双程转换为距离
            detect_amp = envelope(detect_idx);
            
            % 转换为笛卡尔坐标
            % x: 沿航向, y: 跨航向, z: 深度
            x_coord = detect_range * cos(azimuth_rad) * cos(pitch_rad);
            y_coord = detect_range * sin(azimuth_rad);
            z_coord = detect_range * cos(azimuth_rad) * sin(pitch_rad);
            
            % 存储检测结果
            detection_count = detection_count + 1;
            sounding_data.pitch(detection_count) = pitch_angles(p);
            sounding_data.azimuth(detection_count) = azimuth_angles(a);
            sounding_data.range(detection_count) = detect_range;
            sounding_data.amplitude(detection_count) = detect_amp;
            sounding_data.x(detection_count) = x_coord;
            sounding_data.y(detection_count) = y_coord;
            sounding_data.z(detection_count) = z_coord;
        end
    end
end

close(wb);

disp(['海底检测完成! 共检测到 ', num2str(detection_count), ' 个测深点']);

%% 显示波束形成结果 - 单个俯仰角的方位-距离图

disp(' ');
disp('正在生成可视化结果...');

figure('Name', 'MBES 波束形成结果', 'Position', [50 60 1400 600]);

% 选择中间俯仰角进行显示
mid_pitch_idx = round(N_pitch/2);

% 子图1: 时域波束数据
subplot(1,2,1);
range_axis = time_axis * c0 / 2;
beam_image = squeeze(beam_data(mid_pitch_idx, :, :))';
imagesc(azimuth_angles, range_axis, 20*log10(abs(beam_image) + eps));
xlabel('方位角 [度]'); ylabel('距离 [m]');
title(['俯仰角 ', num2str(pitch_angles(mid_pitch_idx)), '° 的方位-距离图']);
colormap(hot); colorbar;
caxis([-60 0]);
axis tight;

% 子图2: 所有俯仰角的最大投影
subplot(1,2,2);
max_projection = squeeze(max(abs(beam_data), [], 3));
imagesc(azimuth_angles, pitch_angles, 20*log10(max_projection + eps));
xlabel('方位角 [度]'); ylabel('俯仰角 [度]');
title('所有波束的最大强度投影');
colormap(hot); colorbar;
caxis([-40 0]);
axis tight;

%% 显示测深点云
% 
% 以3D点云形式显示检测到的海底点

figure('Name', 'MBES 测深点云', 'Position', [100 100 1400 700]);

% 子图1: 3D点云视图
subplot(2,2,1);
scatter3(sounding_data.x, sounding_data.y, sounding_data.z, 20, ...
         sounding_data.amplitude, 'filled');
hold on;
plot3(0, 0, 0, 'rp', 'MarkerFaceColor', 'r', 'MarkerSize', 15);
colormap(jet); colorbar;
set(gca, 'ZDir', 'reverse');
axis equal tight; grid on;
xlabel('沿航向 x [m]'); ylabel('跨航向 y [m]'); zlabel('深度 z [m]');
title('3D 测深点云');
view(45, 30);

% 子图2: 顶视图
subplot(2,2,2);
scatter(sounding_data.x, sounding_data.y, 20, sounding_data.z, 'filled');
hold on;
plot(0, 0, 'rp', 'MarkerFaceColor', 'r', 'MarkerSize', 12);
colormap(jet); colorbar;
axis equal tight; grid on;
xlabel('沿航向 x [m]'); ylabel('跨航向 y [m]');
title('顶视图 (颜色表示深度)');

% 子图3: 侧视图 (沿航向)
subplot(2,2,3);
scatter(sounding_data.x, sounding_data.z, 20, sounding_data.amplitude, 'filled');
hold on;
plot(0, 0, 'rp', 'MarkerFaceColor', 'r', 'MarkerSize', 12);
colormap(jet); colorbar;
set(gca, 'YDir', 'reverse');
axis tight; grid on;
xlabel('沿航向 x [m]'); ylabel('深度 z [m]');
title('沿航向侧视图');

% 子图4: 侧视图 (跨航向)
subplot(2,2,4);
scatter(sounding_data.y, sounding_data.z, 20, sounding_data.amplitude, 'filled');
hold on;
plot(0, 0, 'rp', 'MarkerFaceColor', 'r', 'MarkerSize', 12);
colormap(jet); colorbar;
set(gca, 'YDir', 'reverse');
axis tight; grid on;
xlabel('跨航向 y [m]'); ylabel('深度 z [m]');
title('跨航向侧视图');

%% 生成规则网格地形图
% 
% 将离散测深点插值到规则网格,生成连续海底地形

figure('Name', 'MBES 海底地形图', 'Position', [150 150 1200 800]);

% 创建插值网格
x_grid = linspace(min(sounding_data.x), max(sounding_data.x), 200);
y_grid = linspace(min(sounding_data.y), max(sounding_data.y), 200);
[X_grid, Y_grid] = meshgrid(x_grid, y_grid);

% 网格化深度数据
Z_grid = griddata(sounding_data.x, sounding_data.y, sounding_data.z, ...
                  X_grid, Y_grid, 'natural');

% 子图1: 3D地形
subplot(2,2,1);
surf(X_grid, Y_grid, Z_grid);
hold on;
plot3(0, 0, min(Z_grid(:)), 'rp', 'MarkerFaceColor', 'r', 'MarkerSize', 15);
shading interp;
colormap(jet); colorbar;
set(gca, 'ZDir', 'reverse');
axis tight; grid on;
xlabel('沿航向 x [m]'); ylabel('跨航向 y [m]'); zlabel('深度 z [m]');
title('插值后的3D海底地形');
view(45, 30);

% 子图2: 等深线图
subplot(2,2,2);
contourf(X_grid, Y_grid, Z_grid, 20);
hold on;
plot(0, 0, 'rp', 'MarkerFaceColor', 'r', 'MarkerSize', 12);
colormap(jet); colorbar;
axis equal tight; grid on;
xlabel('沿航向 x [m]'); ylabel('跨航向 y [m]');
title('深度等值线图');

% 子图3: 沿航向中线剖面对比
subplot(2,2,3);
mid_y_idx = round(size(Z_grid, 1)/2);
plot(x_grid, Z_grid(mid_y_idx, :), 'b-', 'LineWidth', 2, 'DisplayName', '测量地形');
hold on;
% 如果有真实地形,绘制对比
if exist('Zg', 'var')
    true_profile = Zg(round(Ny/2), :);
    plot(x_lin, true_profile, 'r--', 'LineWidth', 1.5, 'DisplayName', '真实地形');
end
set(gca, 'YDir', 'reverse');
grid on; legend('show');
xlabel('沿航向 x [m]'); ylabel('深度 z [m]');
title('沿航向中线剖面');

% 子图4: 跨航向中线剖面对比
subplot(2,2,4);
mid_x_idx = round(size(Z_grid, 2)/2);
plot(y_grid, Z_grid(:, mid_x_idx), 'b-', 'LineWidth', 2, 'DisplayName', '测量地形');
hold on;
if exist('Zg', 'var')
    true_profile = Zg(:, round(Nx/2));
    plot(y_lin, true_profile, 'r--', 'LineWidth', 1.5, 'DisplayName', '真实地形');
end
set(gca, 'YDir', 'reverse');
grid on; legend('show');
xlabel('跨航向 y [m]'); ylabel('深度 z [m]');
title('跨航向中线剖面');

%% 统计信息

disp(' ');
disp('========================================');
disp('测深统计信息');
disp('========================================');
disp(['检测点数: ', num2str(detection_count)]);
disp(['覆盖范围 (沿航向): ', num2str(min(sounding_data.x), '%.1f'), ...
      ' - ', num2str(max(sounding_data.x), '%.1f'), ' m']);
disp(['覆盖范围 (跨航向): ', num2str(min(sounding_data.y), '%.1f'), ...
      ' - ', num2str(max(sounding_data.y), '%.1f'), ' m']);
disp(['深度范围: ', num2str(min(sounding_data.z), '%.1f'), ...
      ' - ', num2str(max(sounding_data.z), '%.1f'), ' m']);
disp(['平均深度: ', num2str(mean(sounding_data.z), '%.1f'), ' m']);
disp(['深度标准差: ', num2str(std(sounding_data.z), '%.2f'), ' m']);
disp(' ');
disp('多波束测深声纳仿真完成!');

%% 保存数据 (可选)
% 
% 可以将测深数据保存为MAT文件或导出为其他格式

save_data = false;  % 设置为true以保存数据

if save_data
    disp('正在保存数据...');
    
    % 保存通道数据(UFF格式)
    % channel_data.write('MBES_channel_data.uff', 'channel_data');
    
    % 保存测深点云数据
    save('MBES_sounding_data.mat', 'sounding_data', 'beam_data', ...
         'pitch_angles', 'azimuth_angles', 'X_grid', 'Y_grid', 'Z_grid');
    
    disp('数据保存完成!');
end

%% 导出测深数据为XYZ格式 (可选)
% 
% 导出为通用的点云格式,便于在其他软件中处理

export_xyz = false;  % 设置为true以导出

if export_xyz
    filename = 'MBES_sounding_points.xyz';
    fid = fopen(filename, 'w');
    fprintf(fid, 'X Y Z Amplitude\n');
    for i = 1:length(sounding_data.x)
        fprintf(fid, '%.3f %.3f %.3f %.6f\n', ...
                sounding_data.x(i), sounding_data.y(i), ...
                sounding_data.z(i), sounding_data.amplitude(i));
    end
    fclose(fid);
    disp(['点云数据已导出到: ', filename]);
end