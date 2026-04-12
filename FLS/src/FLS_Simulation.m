%% =========================================================================
% 文件名:    FLS_Simulation.m
% 版本:      2.1 (Path-Optimized)
% 描述:      前视声呐 (FLS) 仿真与波束形成算法对比验证平台
%            已适配模块化目录结构，支持自动路径加载与同目录资源读取
% 作者:      [楼浩/OIT]
% =========================================================================

%% 环境初始化
clear; close all; clc;

% --- 1. 自动处理项目路径 (依赖根目录的 startup.m) ---
% 获取当前脚本所在的绝对路径
current_script_path = fileparts(mfilename('fullpath')); 
% 向上推两级找到项目根目录 (Matlab_code/)
project_root = fullfile(current_script_path, '..', '..'); 

if exist(fullfile(project_root, 'startup.m'), 'file')
    run(fullfile(project_root, 'startup.m'));
else
    % 如果没有 startup.m，尝试手动添加必要的公共工具箱路径
    addpath(genpath(fullfile(project_root, 'common_utils')));
    warning('未找到 startup.m，已尝试手动加载 common_utils 路径。');
end

%% 基本常量定义
c0 = 1500;      % 声速 [m/s]
fs = 10e6;      % 采样频率 [Hz]
dt = 1/fs;      

%% Field II 初始化
field_init(0);
set_field('c', c0);
set_field('fs', fs);
set_field('use_rectangles', 1);

%% 换能器定义 - 128阵元线阵
% 
% 定义前视声纳使用的超声换能器阵列。
% 本示例使用128阵元线阵,阵列沿x轴放置,z为前视方向。

probe = uff.linear_array();
f0                      = 0.9e+06;         % 换能器中心频率 [Hz]
lambda                  = c0/f0;           % 波长 [m]
probe.element_height    = 15e-3;           % 阵元高度 [m]
probe.pitch             = lambda;          % 阵元间距 [m]
kerf                    = 0.1*lambda;      % 阵元间隙 [m]
probe.element_width     = probe.pitch-kerf;% 阵元宽度 [m]
probe.N                 = 128;             % 阵元数量
pulse_duration          = 2.5;             % 脉冲持续时间 [周期]

%% 场景参数定义
% 目标: 使用多个紧邻点目标以展示 CBF 与 FISTA 分辨率对比
% CBF 方位向 3dB 波束宽度 (Blackman 加窗) ≈ 1.68×0.886×λ/D ≈ 0.666°
%   对应 z=5m 处横向分辨率 ≈ 5.8 cm
% 目标间距 4 cm < 5.8 cm → CBF 将三目标融合成一个宽斑, FISTA 可分辨各点

% 海底地形参数
seabed_len      = 15;          % 海底长度(z方向) [m]
seabed_wid      = 15;          % 海底宽度(x方向) [m]
seabed_y0       = 10;          % 海底平均深度 [m]
seabed_amp      = 3;           % 海底起伏幅度 [m]
N_seabed_scatter = 1e4;        % 海底散射点数量

% 视场(FOV)参数
FOV_yaw         = -45:1:45;    % 水平角度网格 [度]
FOV_pitch       = -45:1:45;    % 俯仰角度网格 [度]
FOV_R           = seabed_len;  % 视场射线长度 [m]

%% 脉冲定义
% 
% 使用chirp(线性调频)信号定义脉冲-回波信号。
% 这比简单的高斯脉冲提供更好的分辨率和信噪比。

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

% 显示脉冲以检查延迟估计是否正确(以及脉冲是否对称)
figure;
plot((0:(length(two_way_ir)-1))*dt - lag*dt, two_way_ir); 
hold on; grid on; axis tight
plot((0:(length(two_way_ir)-1))*dt - lag*dt, abs(hilbert(two_way_ir)), 'r')
plot([0 0], [min(two_way_ir) max(two_way_ir)], 'g');
legend('双程脉冲', '包络', '估计延迟');
title('Field II 双程冲激响应');

%% 孔径对象
% 
% 使用Field II的*xdc_linear_array*函数定义网格几何。

noSubAz = round(probe.element_width/(lambda/8));   % 方位向子阵元数
noSubEl = round(probe.element_height/(lambda/8));  % 俯仰向子阵元数
Th = xdc_linear_array(probe.N, probe.element_width, probe.element_height, kerf, noSubAz, noSubEl, [0 0 Inf]); 
Rh = xdc_linear_array(probe.N, probe.element_width, probe.element_height, kerf, noSubAz, noSubEl, [0 0 Inf]); 

% 设置激励、冲激响应和障板:
xdc_excitation(Th, excitation);
xdc_impulse(Th, impulse_response);
xdc_baffle(Th, 0);
xdc_center_focus(Th, [0 0 0]);
xdc_impulse(Rh, impulse_response);
xdc_baffle(Rh, 0);
xdc_center_focus(Rh, [0 0 0]);

%% 目标点定义
% 坐标格式: [x(横向/方位), y(俯仰/仰角=0 表示在扫描平面内), z(距离)] —— Field II 约定
%
% 3 个紧邻点目标, 间距 4 cm, 全部位于 z = 5 m 处:
%   相邻目标角度间距 ≈ atan(0.04/5) = 0.46°  <  CBF 3dB 宽度 0.67°
%   → CBF: 三目标融合为单个宽斑
%   → FISTA: 去卷积后可分辨三个独立峰值

pos_target = [
    -0.04,  0,  5;   % 目标 1: 左偏 4 cm, 距离 5 m
     0.00,  0,  5;   % 目标 2: 正中,      距离 5 m
     0.04,  0,  5;   % 目标 3: 右偏 4 cm, 距离 5 m
];

disp(['目标点定义完成: ', num2str(size(pos_target,1)), ' 个紧邻散射点 (间距 4 cm @ 5 m)']);

%% 定义海底散射体
% 
% 生成具有起伏的海底地形散射点

disp('正在生成海底散射体...');

% 计算网格点数
Nx = round(sqrt(N_seabed_scatter * seabed_wid / seabed_len));
Nz = round(N_seabed_scatter / Nx);

% 生成均匀网格
x_lin = linspace(-seabed_wid/2, seabed_wid/2, Nx);
z_lin = linspace(0, seabed_len, Nz);
[Xg, Zg] = meshgrid(x_lin, z_lin);

% 添加深度起伏(y为深度方向)
Yg = seabed_y0 + seabed_amp * (0.5*sin(2*pi*Xg/seabed_wid) + ...  % 横向波动
                                0.5*sin(2*pi*Zg/seabed_len)) + ...  % 纵向波动
                                0.1*randn(size(Xg));                % 随机噪声


% 展开为点列表
pos_seabed = [Xg(:) Yg(:) Zg(:)];

% 设置散射系数
amp_target = 0.80 * ones(size(pos_target,1), 1);      % 目标强反射 (高于海底, 便于对比)
amp_seabed = 0.10 * abs(randn(size(pos_seabed,1), 1)); % 海底弱随机散射

disp(['海底散射体生成完成: ', num2str(size(pos_seabed,1)), ' 个散射点']);

%% 场景可视化
% 
% 显示完整的3D场景,包括目标、海底和视场网格

disp('正在生成场景预览...');

figure('Name', '前视声纳场景预览', 'Position', [50 60 1400 700]);

% 子图1: 3D视图
subplot(2,3,1);
scatter3(pos_seabed(:,1), pos_seabed(:,2), pos_seabed(:,3), 3, 'b.'); 
hold on;
scatter3(pos_target(:,1), pos_target(:,2), pos_target(:,3), 10, 'g', 'filled');
plot3(0, 0, 0, 'rp', 'MarkerFaceColor', 'r', 'MarkerSize', 10);
set(gca, 'ZDir', 'reverse');
axis equal tight; grid on;
xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
title('3D视图');
view(35, 25);

% 子图2: 顶视图(x-y)
subplot(2,3,2);
scatter(pos_seabed(:,1), pos_seabed(:,2), 3, 'b.'); 
hold on;
scatter(pos_target(:,1), pos_target(:,2), 10, 'g', 'filled');
plot(0, 0, 'rp', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
axis equal tight; grid on;
xlabel('x [m]'); ylabel('y [m]');
title('顶视图 (x-y)');

% 子图3: 侧视图(y-z)
subplot(2,3,3);
scatter(pos_seabed(:,2), pos_seabed(:,3), 3, 'b.'); 
hold on;
scatter(pos_target(:,2), pos_target(:,3), 10, 'g', 'filled');
plot(0, 0, 'rp', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
set(gca, 'YDir', 'reverse');
axis equal tight; grid on;
xlabel('y [m]'); ylabel('z [m]');
title('侧视图 (y-z)');

% 子图4-6: 视场网格和散射体
subplot(2,3,[4 5 6]);
hold on; axis equal; grid on;
set(gca, 'ZDir', 'reverse');
xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
title('视场网格与散射体');

% 绘制水平角度网格线
for ya = FOV_yaw
    dir = [sind(ya) cosd(ya) 0] * FOV_R;
    plot3([0 dir(1)], [0 0], [0 dir(2)], 'b-');
end

% 绘制俯仰角度网格线
for pa = FOV_pitch
    dir = [0 cosd(pa) sind(pa)] * FOV_R;
    plot3([0 dir(1)], [0 dir(3)], [0 dir(2)], 'r-');
end

% 绘制稀疏采样的散射体(提高显示速度)
scatter3(pos_seabed(1:15:end,1), pos_seabed(1:15:end,2), ...
         pos_seabed(1:15:end,3), 3, 'b.');
scatter3(pos_target(:,1), pos_target(:,2), pos_target(:,3), 10, 'g', 'filled');
plot3(0, 0, 0, 'rp', 'MarkerFaceColor', 'r', 'MarkerSize', 10);
legend('水平角网格', '俯仰角网格', '海底散射体', '目标', '阵列中心', 'Location', 'best');
view(3);
drawnow;



%% 定义发射序列
% 
% 定义单个平面波发射(无角度扫描)

F_number = 1;           % F数 = 焦距/孔径宽度
alpha_max = 0;          % 最大发射角度 [rad]
Na = 1;                 % 平面波数量
F = 1;                  % 帧数
alpha = linspace(-alpha_max, alpha_max, Na);  % 角度向量 [rad]

%% 组合所有散射点
% 
% 将目标和海底散射点合并为完整的散射场

point_position = [pos_target; pos_seabed];
point_amplitudes = [amp_target; amp_seabed];

disp(['总散射点数: ', num2str(size(point_position,1))]);


%% ========================================================================
%  新增模块：绘制前视声纳(FLS) 3D 探测范围 (符合 FLS 视觉习惯版)
%  坐标映射说明：
%  - 绘图 X轴 = 物理横向 (Field II X)
%  - 绘图 Y轴 = 物理前方距离 (Field II Z)  <-- 关键交换
%  - 绘图 Z轴 = 物理深度 (Field II Y)      <-- 关键交换
% ========================================================================
figure(99); clf;
set(gcf, 'Name', 'FLS 3D场景与FOV可视化', 'Color', 'w');
hold on; box on; grid on;

% --- 1. 获取探测几何参数 ---
R_max_vis = seabed_len; % 最大探测距离 (15m)
yaw_half_deg = max(abs(FOV_yaw));   
pitch_half_deg = max(abs(FOV_pitch)); 

% --- 2. 转换点云坐标 (关键步骤：交换 Y 和 Z) ---
% point_position 原本是 [x, y(深度), z(距离)]
% 我们将其转换为绘图用的 [x, z(距离), y(深度)]
plot_pos_scat = point_position; 
plot_pos_scat(:,2) = point_position(:,3); % 绘图Y = 物理Z
plot_pos_scat(:,3) = point_position(:,2); % 绘图Z = 物理Y

% --- 3. 计算视锥体顶点 (已适配交换后的坐标系) ---
% Field II 物理尺寸
X_face_half = R_max_vis * tand(yaw_half_deg);
Y_face_half = R_max_vis * tand(pitch_half_deg); % 这是物理上的半高

% 顶点定义 [绘图X, 绘图Y(距离), 绘图Z(深度)]
V0 = [0, 0, 0]; % 原点
% 最远端面的四个角点 (距离都在 R_max_vis)
V1 = [-X_face_half, R_max_vis, -Y_face_half];
V2 = [ X_face_half, R_max_vis, -Y_face_half];
V3 = [ X_face_half, R_max_vis,  Y_face_half];
V4 = [-X_face_half, R_max_vis,  Y_face_half];

Vertices_pyramid = [V0; V1; V2; V3; V4];

% --- 4. 绘制半透明视锥体 ---
Faces_pyramid = [
    1 2 3 NaN;  % 侧面
    1 3 4 NaN;  % 侧面
    1 4 5 NaN;  % 侧面
    1 5 2 NaN;  % 侧面
    2 3 4 5     % 底面 (最远端)
];

patch('Vertices', Vertices_pyramid, 'Faces', Faces_pyramid, ...
      'FaceColor', 'g', 'FaceAlpha', 0.1, ...
      'EdgeColor', 'g', 'LineWidth', 1.5, 'LineStyle', '--', ...
      'DisplayName', 'FOV探测边界');

% --- 5. 绘制场景物体 ---
% 声纳头
plot3(0, 0, 0, 'ks', 'MarkerSize', 12, 'MarkerFaceColor', 'k', 'DisplayName', '声纳头');

% 场景点云 (使用转换后的坐标 plot_pos_scat)
scatter3(plot_pos_scat(:,1), plot_pos_scat(:,2), plot_pos_scat(:,3), ...
         12, point_amplitudes, 'filled', 'DisplayName', '海底与目标');
colormap(jet); 

% --- 6. 视觉调整 (符合航海习惯) ---
xlabel('横向 X [m]');
ylabel('前方距离 Z [m]');  % 注意：现在Y轴显示的是距离
zlabel('深度 Y [m]');      % 注意：现在Z轴显示的是深度

% 【核心设置】反转 Z 轴
% 因为深度通常是向下增加的，但在 3D 图里 Z 轴向上是正。
% 我们反转它，让 Z=0 (水面) 在最上面，Z=10 (海底) 在下面。
set(gca, 'ZDir', 'reverse'); 

% 调整视角
% view(-30, 30); % 侧俯视
view(0, 90); % 切换到【顶视图/地图视角】 (X-Y平面)，就像导航地图一样
% 如果想看立体效果，可以取消上面一行，改用: view(45, 30);

axis equal tight;
title('FLS 仿真场景 (符合直觉坐标系)');
legend('show', 'Location', 'northeast');

fprintf('已更新 FLS 可视化窗口 (Figure 99)，坐标轴已修正。\n');
%% 输出数据初始化

cropat = round(2*50e-3/c0/dt);    % 最大时间采样点
CPW = zeros(cropat, probe.N, Na, F);  % 通道数据矩阵

%% 计算通道数据
% 
% 使用Field II计算每个发射角度的散射响应

time_index = 0;
disp('Field II: 正在计算前视声纳数据集');
wb = waitbar(0, 'Field II: 正在计算前视声纳数据集');

% Blackman窗用于接收孔径加权
win = blackman(probe.N);

for f = 1:F
    waitbar(0, wb, sprintf('Field II: 正在计算数据集, 第 %d 帧', f));
    
    for n = 1:Na
        waitbar(n/Na, wb);
        disp(['正在计算角度 ', num2str(n), ' / ', num2str(Na)]);
        
        % 发射孔径设置
        tmp = zeros(1, probe.N/2);
        xdc_apodization(Th, 0, [tmp(1:end-1), 1, 1, tmp(1:end-1)]);
        xdc_times_focus(Th, 0, probe.geometry(:,1)' .* sin(alpha(n))/c0);
        
        % 接收孔径设置
        xdc_apodization(Rh, 0, ones(1, probe.N));
        xdc_focus_times(Rh, 0, zeros(1, probe.N));
        
        % 执行散射计算
        [v, t] = calc_scat_multi(Th, Rh, point_position, point_amplitudes);
        
        % 匹配滤波 + 加窗
        % 用于降低旁瓣,改善图像对比度
        sig_mf = zeros(size(v));
        for m = 1:probe.N
            tmp_sig = conv(v(:,m), MF, 'same');
            sig_mf(:,m) = tmp_sig * win(m);
        end
        
        % 存储数据
        CPW(1:size(v,1), :, n, f) = sig_mf;
        
        % 保存发射序列信息
        seq(n) = uff.wave();
        seq(n).probe = probe;
        seq(n).source.distance = Inf;  % 平面波
        seq(n).sound_speed = c0;
        seq(n).delay = -lag*dt + t;
    end
end
close(wb);

disp('数据计算完成!');

%% 通道数据
% 
% 创建UFF数据结构来存储捕获的超声通道数据

channel_data = uff.channel_data();
channel_data.sampling_frequency = fs;
channel_data.sound_speed = c0;
channel_data.initial_time = 0;
channel_data.pulse = pulse;
channel_data.probe = probe;
channel_data.sequence = seq;
channel_data.data = CPW ./ max(CPW(:));  % 归一化

%% 扫描区域定义
%
% 扫描区域定义为覆盖感兴趣区域的像素集合。
% 这里使用*sector_scan*结构,用方位角和深度定义。

depth_axis = linspace(0, 25, 512)';           % 深度轴: 0-25m, 512点
azimuth_axis = linspace(-45, 45, 451)/180*pi; % 方位轴: ±45°, 451点

sca = uff.sector_scan('azimuth_axis', azimuth_axis, 'depth_axis', depth_axis);

%% 波束形成管道
%
% 有了*channel_data*和*scan*,我们就有了生成超声图像所需的一切。
% 现在使用USTB的*pipeline*结构,它需要*apodization*结构
% 以及*channel_data*和*scan*。

%% ========================================================================
%  多算法波束形成 (CBF & MVDR) - 最终修正版
% ========================================================================

% 初始化管道
pipe = pipeline();
pipe.channel_data = channel_data;
pipe.scan = sca;
% 设置加窗
pipe.receive_apodization.window = uff.window.tukey25;
pipe.receive_apodization.f_number = F_number;

% -------------------------------------------------------------------------
% 1. 常规波束形成 (CBF / DAS) - 作为基准
% -------------------------------------------------------------------------
disp('正在执行 CBF (DAS) 波束形成...');
% 使用默认 DAS，它会自动进行延时和全孔径求和
proc_cbf = midprocess.das();
b_data_cbf = pipe.go({proc_cbf});

% -------------------------------------------------------------------------
% 2. 自适应波束形成 (MVDR / Capon)
% -------------------------------------------------------------------------
disp('正在执行 MVDR (Capon) 预处理...');

% [关键步骤 A] 预波束形成 (Pre-beamforming)
% 我们需要先用 DAS 计算延时(Focusing)，但保留接收通道的数据不求和。
% 这样 Capon 才能拿到每个像素点对应的、已对齐的通道数据来计算权重。
proc_pre = midprocess.das();
% 设置维度为 dimension.transmit (仅发射求和)，这意味着接收通道(Receive)会被保留
proc_pre.dimension = dimension.transmit; 

% 执行预处理，得到中间数据 (Pixels x Channels)
b_data_pre = pipe.go({proc_pre});

disp('正在执行 MVDR (Capon) 核心计算...');
disp('注意：计算涉及矩阵求逆，速度较慢，请耐心等待...');

% [关键步骤 B] 调用 Capon 后处理
mvdr = postprocess.capon_minimum_variance();
mvdr.dimension = dimension.receive;  % 在接收维上进行自适应加权
mvdr.channel_data = channel_data;    % 提供探头几何信息
mvdr.scan = sca;

% 设置 MVDR 参数
% L_elements: 子阵列长度 (Subarray size)。通常取阵元数的一半左右。
% 这里的 probe.N 是 128，我们设为 64。较小的值计算快但分辨率低，较大的值分辨率高但更慢且不稳定。
mvdr.L_elements = channel_data.probe.N / 2; 

% K_in_lambda: 时间平均平滑因子 (Temporal averaging)
mvdr.K_in_lambda = 1.5; 

% regCoef: 对角加载系数 (正则化)，防止矩阵不可逆
mvdr.regCoef = 1e-2;

% [核心修正] 将预处理后的"聚焦数据"喂给 MVDR
mvdr.input = b_data_pre; 

% 执行 MVDR
b_data_mvdr = mvdr.go();


%% ========================================================================
%  CBF + FISTA 去卷积波束形成
%  算法: Fast Iterative Shrinkage-Thresholding Algorithm (FISTA)
%  参考: Beck & Teboulle (2009); Sonar sparse deconvolution beamforming
%
%  模型: y = H(x) + n
%    y: CBF 包络图像 (观测值)
%    H: 点扩散函数算子 (PSF 循环卷积)
%    x: 待恢复的稀疏反射率图
%    n: 噪声
%
%  目标函数: min_{x >= 0}  (1/2)*||H(x) - y||_F^2 + lambda*||x||_1
% ========================================================================

disp('正在执行 CBF + FISTA 去卷积...');

%% Step 1: 从 CBF 波束形成数据中提取 2D 复数图像
N_d = numel(depth_axis);    % 深度像素数: 512
N_a = numel(azimuth_axis);  % 方位像素数: 451

% b_data_cbf.data 存储为 [N_pixels × N_tx], N_pixels = N_d * N_a
% reshape 为 [深度 × 方位] 二维矩阵
raw_cbf = b_data_cbf.data;
data_sz = size(raw_cbf);
if data_sz(1) == N_d && data_sz(2) == N_a
    % 已是 2D 矩阵格式
    cbf_complex = raw_cbf;
elseif data_sz(1) == N_d * N_a
    % 展平格式: [N_pixels × N_tx]
    cbf_complex = reshape(raw_cbf(:, 1), N_d, N_a);
else
    % 三维格式: [N_d × N_a × N_tx]
    cbf_complex = raw_cbf(:, :, 1);
end

cbf_env = abs(cbf_complex);                     % 包络检波
cbf_img = cbf_env / (max(cbf_env(:)) + eps);    % 归一化至 [0, 1]

%% Step 2: 估计系统点扩散函数 (PSF)
% PSF 由阵列方向图 (方位向) 和脉冲带宽 (距离向) 共同决定

% --- 方位向 PSF: Blackman 加窗的阵列方向图 ---
% 均匀线阵 3dB 波束宽度 [rad]: 0.886 * lambda / D_array
bw_az_rad   = 0.886 * lambda / (probe.N * probe.pitch);
bw_az_rad   = bw_az_rad * 1.68;   % Blackman 窗主瓣展宽因子 (~1.68×)
az_step_rad = (azimuth_axis(end) - azimuth_axis(1)) / (N_a - 1);
% 将 FWHM 转为像素域高斯标准差
sigma_az_px = bw_az_rad / az_step_rad / (2 * sqrt(2 * log(2)));

% --- 距离向 PSF: 匹配滤波后的脉冲压缩响应 ---
range_res_m  = c0 / (2 * BW);     % 距离分辨率 [m]
depth_step_m = (depth_axis(end) - depth_axis(1)) / (N_d - 1);
sigma_r_px   = range_res_m / depth_step_m / (2 * sqrt(2 * log(2)));
sigma_r_px   = max(sigma_r_px, 0.5);  % 至少保留 0.5 像素宽

fprintf('  PSF: σ_方位 = %.2f px (%.2f°),  σ_距离 = %.2f px (%.1f mm)\n', ...
    sigma_az_px, sigma_az_px * (2*sqrt(2*log(2))) * az_step_rad * 180/pi, ...
    sigma_r_px,  sigma_r_px  * (2*sqrt(2*log(2))) * depth_step_m * 1e3);

% --- 构建 2D 高斯 PSF 核 ---
half_w = max(ceil(4 * max(sigma_az_px, sigma_r_px)), 8);
[gx, gy] = meshgrid(-half_w:half_w, -half_w:half_w);
% gx: 方位向 (列方向), gy: 深度向 (行方向)
PSF = exp(-0.5 * ((gx / sigma_az_px).^2 + (gy / sigma_r_px).^2));
PSF = PSF / sum(PSF(:));    % 总能量归一化

%% Step 3: 零填充 (避免循环卷积绕回伪影)
% 参考: EMBED (dtu-act/EMBED) zeropad 策略
% 原始循环卷积将图像边界与对侧"接续", 产生伪影.
% 在图像四周各填充 PSF 半宽 (pad_r / pad_a) 个零像素后再做 FISTA,
% 收敛后裁剪回原始尺寸即得到无绕回伪影的去卷积结果.

Np_r  = size(PSF, 1);          % PSF 深度方向像素数
Np_a  = size(PSF, 2);          % PSF 方位方向像素数
pad_r = floor(Np_r / 2);       % 深度方向单侧填充量
pad_a = floor(Np_a / 2);       % 方位方向单侧填充量
N_pad_d = N_d + 2 * pad_r;     % 填充后深度尺寸
N_pad_a = N_a + 2 * pad_a;     % 填充后方位尺寸

% 观测图像: 居中放入零填充矩阵
cbf_pad = zeros(N_pad_d, N_pad_a);
cbf_pad(pad_r+1 : pad_r+N_d, pad_a+1 : pad_a+N_a) = cbf_img;

% PSF: 放入填充矩阵后 ifftshift → 将 PSF 中心移至 (1,1) (FFT 原点),
% 保证 ifft2 后卷积结果与图像对齐, 无需在迭代内额外调用 fftshift.
PSF_pad = zeros(N_pad_d, N_pad_a);
PSF_pad(1:Np_r, 1:Np_a) = PSF;
PSF_pad = ifftshift(PSF_pad);

%% Step 4: 改进的 FISTA 迭代去卷积
% 算法结构参考 EMBED (dtu-act/EMBED) 的 FISTA.m,
% 在原始 Beck & Teboulle (2009) 框架上做如下改进:
%   1. 零填充消除绕回伪影           (见 Step 3)
%   2. 幂迭代精确估计 Lipschitz 常数 (EMBED lipschitz 函数, 10 次迭代)
%   3. 目标函数历史记录              (EMBED info.obj)
%   4. 耗时统计                      (EMBED info.time)

lambda_fista   = 0.01;   % L1 正则化强度 (越大越稀疏)
max_iter_fista = 200;    % 最大迭代次数

% 预计算填充后 PSF 的 FFT 及伴随 FFT
% 对实数对称 PSF: conj(H_fft) 等价于 EMBED 的 fft2(rot90(PSF,2))
H_fft  = fft2(PSF_pad);
Ht_fft = conj(H_fft);

% 幂迭代估计 Lipschitz 常数 L = λ_max(H^T H)
% 参考: EMBED/src/FISTA.m 内部 lipschitz 函数 (10 次幂迭代)
% 零填充后算子不再是纯循环矩阵, max|H_fft|^2 可能偏高或偏低,
% 幂迭代给出更准确的谱范数估计, 保证收敛且不过于保守.
fprintf('  幂迭代估计 Lipschitz 常数 (10 次迭代)...\n');
v = rand(N_pad_d, N_pad_a);
for pi_k = 1:10
    Av   = real(ifft2(H_fft  .* fft2(v)));
    AtAv = real(ifft2(Ht_fft .* fft2(Av)));
    v_norm = norm(AtAv(:));
    if v_norm < eps; break; end
    v = AtAv / v_norm;
end
L_lip = norm(real(ifft2(H_fft .* fft2(v))), 'fro')^2;
L_lip = max(L_lip, eps);   % 防止除零

fprintf('  FISTA 参数: lambda = %.4f, L = %.4e, 最大迭代 = %d\n', ...
        lambda_fista, L_lip, max_iter_fista);

% FISTA 初始化 (以零填充 CBF 图像为热启动)
x_k          = cbf_pad;
z_k          = x_k;
t_k          = 1;
obj_history  = nan(max_iter_fista, 1);   % 目标函数历史 (参考 EMBED info.obj)
prev_obj     = inf;
t_fista_start = tic;

for k = 1:max_iter_fista
    % 1. 前向卷积: H * z_k
    Hz = real(ifft2(H_fft .* fft2(z_k)));

    % 2. 残差: H*z_k - y
    residual = Hz - cbf_pad;

    % 3. 伴随梯度: H^T * residual
    grad = real(ifft2(Ht_fft .* fft2(residual)));

    % 4. 梯度下降步
    u = z_k - grad / L_lip;

    % 5. 近端算子: 软阈值 (L1 prox) + 非负投影
    tau   = lambda_fista / L_lip;
    x_new = max(u - tau, 0);

    % 6. FISTA Nesterov 动量更新 (参考 EMBED FISTA.m)
    t_new = (1 + sqrt(1 + 4 * t_k^2)) / 2;
    z_k   = x_new + ((t_k - 1) / t_new) * (x_new - x_k);

    % 7. 目标函数记录 & 收敛检测 (参考 EMBED info.obj)
    % x >= 0, 故 ||x||_1 = sum(x)
    obj_val        = 0.5 * norm(residual(:))^2 + lambda_fista * sum(x_new(:));
    obj_history(k) = obj_val;
    rel_change     = abs(prev_obj - obj_val) / (abs(prev_obj) + eps);
    if rel_change < 1e-6 && k > 10
        fprintf('  FISTA 收敛于第 %d 次迭代 (Δobj = %.2e)\n', k, rel_change);
        obj_history = obj_history(1:k);
        break;
    end
    if mod(k, 50) == 0
        fprintf('  迭代 %3d/%d  代价函数 = %.6f\n', k, max_iter_fista, obj_val);
    end

    prev_obj = obj_val;
    x_k      = x_new;
    t_k      = t_new;
end

% 耗时统计 (参考 EMBED info.time)
t_fista_elapsed = toc(t_fista_start);
fprintf('  FISTA 耗时: %.2f 秒\n', t_fista_elapsed);

% 从零填充矩阵裁剪回原始图像尺寸
img_fista = x_k(pad_r+1 : pad_r+N_d, pad_a+1 : pad_a+N_a);
disp('FISTA 去卷积完成!');


%% ========================================================================
%  三路算法对比可视化: CBF / MVDR / CBF+FISTA
% ========================================================================

az_deg = azimuth_axis * (180 / pi);   % 弧度 → 度 (用于坐标轴标注)

fig_cmp = figure('Name', '前视声纳成像算法三路对比 (CBF / MVDR / FISTA)', ...
                 'Position', [50, 50, 1800, 600], 'Color', 'w');

% 三目标位于 z=5m, x∈[-4cm, +4cm] → 方位角 ≈ ±0.46°
% 统一坐标轴范围: 方位 ±5°, 深度 3-8 m
% (海底在 10-18 m 处不会出现在此窗口, 画面干净, 便于分辨率对比)
az_lim    = [-5, 5];   % 方位角显示范围 [°]
depth_lim = [3, 8];    % 深度显示范围 [m]

% --- 子图 1: CBF / DAS ---
ax1 = subplot(1, 3, 1);
b_data_cbf.plot(ax1);
colormap(ax1, hot);
clim([-60 0]);
title('CBF / DAS  (常规波束形成)', 'FontSize', 13, 'FontWeight', 'bold');
xlabel('方位角 [°]'); ylabel('深度 [m]');
set(ax1, 'XDir', 'normal', 'XLim', az_lim, 'YLim', depth_lim);

% --- 子图 2: MVDR / Capon ---
ax2 = subplot(1, 3, 2);
b_data_mvdr.plot(ax2);
colormap(ax2, hot);
clim([-60 0]);
title('MVDR / Capon  (自适应波束形成)', 'FontSize', 13, 'FontWeight', 'bold');
xlabel('方位角 [°]'); ylabel('深度 [m]');
set(ax2, 'XDir', 'normal', 'XLim', az_lim, 'YLim', depth_lim);

% --- 子图 3: CBF + FISTA 去卷积 ---
ax3 = subplot(1, 3, 3);
fista_db = 20 * log10(img_fista / (max(img_fista(:)) + eps) + eps);
fista_db = max(fista_db, -60);        % 限制动态范围到 -60 dB
imagesc(az_deg, depth_axis, fista_db);
colormap(ax3, hot);
clim([-60 0]);
colorbar;
title('CBF + FISTA 去卷积  (稀疏正则化)', 'FontSize', 13, 'FontWeight', 'bold');
xlabel('方位角 [°]'); ylabel('深度 [m]');
set(ax3, 'XDir', 'normal', 'XLim', az_lim, 'YLim', depth_lim);

% 联动三轴: 交互缩放/平移时同步
linkaxes([ax1, ax2, ax3], 'xy');

sgtitle(sprintf('前视声纳成像算法对比  —  3 紧邻目标 (间距 4 cm @ 5 m, CBF 分辨率 %.1f cm)\n CBF 无法分辨 | MVDR 边缘改善 | FISTA 完全分辨', ...
        5 * tand(0.886 * (180/pi) * lambda / (probe.N * probe.pitch) * 1.68) * 100), ...
        'FontSize', 13, 'FontWeight', 'bold');
drawnow;

fprintf('\n=== 前视声纳三路算法对比完成 ===\n');
fprintf('  %-12s 旁瓣高, 分辨率受限于阵列孔径\n', 'CBF / DAS:');
fprintf('  %-12s 自适应抑制旁瓣, 近场效果好, 计算量大\n', 'MVDR:');
fprintf('  %-12s 稀疏去卷积, 主瓣最窄, 目标定位精度最高\n', 'FISTA:');

disp('算法对比完成！');

% 仿真结束，释放 Field II 资源
field_end();

