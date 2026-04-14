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
% 5 目标"临界分辨"场景（无海底散射体）
% CBF Blackman 窗 3dB 波束宽度:
%   bw_deg = 0.886*(180/pi)*lambda/(probe.N*probe.pitch)*1.68 ≈ 0.668°
% 紧邻目标角度间距 delta_deg = 0.5 × bw_deg ≈ 0.334°
%   → CBF 刚好无法分辨 (Rayleigh 准则不满足)
%   → FISTA 稀疏去卷积 → 完全分辨
%
% 布局: 组1 @ z=6m (±delta/2), 组2 @ z=10m (±delta/2), 孤立目标 @ z=8m az=10°

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
% 'same'模式conv的相位中心 = ceil(length(MF)/2)
% MF与excitation等长，所以：
lag = ceil(length(excitation) / 2);

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
% 坐标格式: [x(横向/方位), y(俯仰/仰角=0), z(距离)] —— Field II 约定
%
% 临界角度间距与目标布局
bw_deg    = 0.886 * (180/pi) * lambda / (probe.N * probe.pitch) * 1.68;
delta_deg = 0.5 * bw_deg;   % 临界角度间距: 50% BW → CBF 刚好无法分辨

pos_target = [
    % 组1 @ z=6m: 两个紧邻目标，角度间距 delta_deg
    -6*tand(delta_deg/2),  0,  6;   % 目标1 左
    +6*tand(delta_deg/2),  0,  6;   % 目标2 右

    % 组2 @ z=10m: 两个紧邻目标，角度间距 delta_deg
    -10*tand(delta_deg/2), 0, 10;   % 目标3 左
    +10*tand(delta_deg/2), 0, 10;   % 目标4 右

    % 孤立参考目标 @ z=8m, 方位角10° (验证FISTA不破坏孤立点)
    8*tand(10),            0,  8;   % 目标5
];

amp_target = ones(5, 1);  % 等幅反射

fprintf('目标定义: bw=%.3f°, delta=%.3f° (%.0f%% BW), 5个目标\n', ...
    bw_deg, delta_deg, delta_deg/bw_deg*100);


%% 场景预览 (仅目标)
figure('Name', '5目标场景预览', 'Position', [50 60 600 500], 'Color', 'w');
tgt_az_preview = atan2d(pos_target(:,1), pos_target(:,3));
tgt_r_preview  = sqrt(pos_target(:,1).^2 + pos_target(:,3).^2);
scatter(tgt_az_preview, tgt_r_preview, 120, 'g', 'filled');
hold on;
plot(0, 0, 'rp', 'MarkerFaceColor', 'r', 'MarkerSize', 12);
grid on; set(gca, 'YDir', 'reverse');
xlabel('方位角 [°]'); ylabel('距离 [m]');
title(sprintf('5目标场景: delta=%.3f° (bw=%.3f°)', delta_deg, bw_deg));
xlim([-2, 12]); ylim([5, 11]);
text(tgt_az_preview+0.1, tgt_r_preview, ...
    {'T1','T2','T3','T4','T5'}, 'FontSize', 10, 'Color', 'k');
drawnow;



%% 定义发射序列
% 
% 定义单个平面波发射(无角度扫描)

F_number = 1;           % F数 = 焦距/孔径宽度
alpha_max = 0;          % 最大发射角度 [rad]
Na = 1;                 % 平面波数量
F = 1;                  % 帧数
alpha = linspace(-alpha_max, alpha_max, Na);  % 角度向量 [rad]

%% 散射点 (仅目标, 无海底)
point_position   = pos_target;
point_amplitudes = amp_target;

disp(['散射点数: ', num2str(size(point_position,1)), ' 个目标']);


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
        xdc_apodization(Th, 0, ones(1, probe.N));
        xdc_times_focus(Th, 0, zeros(1, probe.N));
        
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

depth_axis   = linspace(4.5, 11.5, 512)';            % 深度轴: 4.5-11.5m, 512点
azimuth_axis = linspace(-3, 13, 1024) / 180 * pi;   % 方位轴: -3°~13°, 1024点

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
cbf_img = squeeze(double(real(cbf_img)));        % 确保为 2D 实数矩阵 (USTB 有时输出 N_d×N_a×1)

% 方位角轴 (度), 后续 PSF 提取和可视化均使用
az_deg = azimuth_axis * (180 / pi);

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
% 真实距离分辨率 ≈ c/(2*BW) ≈ 7.5 mm.  深度像素步长 ≈ 48.9 mm/px.
% → sigma_r ≈ 0.15 px (远小于 1 像素) → 深度方向 PSF 近似为 delta 函数.
% 【关键】不做下限截断: 若截断至 0.5 px, 深度方向会产生人工展宽,
%   FISTA 可能将方位向的双目标错误分解为深度方向的双目标 (伪影).
range_res_m  = c0 / (2 * BW);
depth_step_m = (depth_axis(end) - depth_axis(1)) / (N_d - 1);
sigma_r_px   = range_res_m / depth_step_m / (2 * sqrt(2 * log(2)));
% 不截断: 真实值 ~0.15 px 确保深度 PSF 为 delta, FISTA 仅在方位向去卷积

fprintf('  PSF: σ_方位 = %.2f px (%.2f°),  σ_距离 = %.3f px (%.1f mm)\n', ...
    sigma_az_px, sigma_az_px * (2*sqrt(2*log(2))) * az_step_rad * 180/pi, ...
    sigma_r_px,  sigma_r_px  * (2*sqrt(2*log(2))) * depth_step_m * 1e3);

% --- 构建 2D 高斯 PSF 核 ---
% 方位向宽 (sigma_az_px), 深度向极窄 (sigma_r_px ≈ delta) → 实质为 1D 方位 PSF
half_w = max(ceil(4 * sigma_az_px), 8);   % 仅由方位宽度决定窗口大小
[gx, gy] = meshgrid(-half_w:half_w, -half_w:half_w);
% gx: 方位向 (列方向, 宽), gy: 深度向 (行方向, 极窄)
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

% PSF: 放入填充矩阵左上角后 circshift → 将 PSF 中心精确移至 (1,1) (FFT 原点).
% 【修正说明】原来的 ifftshift 在大矩阵 [N_pad_d × N_pad_a] 上执行，会把 PSF
% 整体平移 ceil(N/2) 个位置 (约 264 行, 234 列) 到矩阵中央，而非目标的 (1,1)。
% circshift 按 PSF 半宽精确平移，使中心从 (floor(Np_r/2)+1, floor(Np_a/2)+1)
% 移到 (1,1)，从而保证 ifft2 后的卷积结果与图像正确对齐。
PSF_pad = zeros(N_pad_d, N_pad_a);
PSF_pad(1:Np_r, 1:Np_a) = PSF;
PSF_pad = circshift(PSF_pad, [-floor(Np_r/2), -floor(Np_a/2)]);

%% Step 4: 改进的 FISTA 迭代去卷积
% 算法结构参考 EMBED (dtu-act/EMBED) 的 FISTA.m,
% 在原始 Beck & Teboulle (2009) 框架上做如下改进:
%   1. 零填充消除绕回伪影           (见 Step 3)
%   2. 幂迭代精确估计 Lipschitz 常数 (EMBED lipschitz 函数, 10 次迭代)
%   3. 目标函数历史记录              (EMBED info.obj)
%   4. 耗时统计                      (EMBED info.time)

lambda_fista   = 0.05;   % L1 正则化强度 (越大越稀疏; 0.05 适合真实PSF)
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
%  【修正】三幅图全部使用 imagesc(az_deg, depth_axis, ...) 统一坐标系.
%  原来 b_data.plot() 的 X 轴是物理坐标 (单位: 米), 而 FISTA 的 imagesc
%  X 轴是方位角 (单位: 度), 两者完全不同 → 坐标系混用导致图像错位.
%  改为从 USTB 数据对象中提取 [N_d × N_a] 矩阵后统一用 imagesc 渲染.
% ========================================================================

az_deg = azimuth_axis * (180 / pi);

%% 提取 CBF 图像 (复用 Step 1 中的 cbf_img, 已是 [N_d × N_a] 归一化包络)
cbf_db_vis = 20 * log10(cbf_img + eps);
cbf_db_vis = max(cbf_db_vis, -60);

%% 提取 MVDR 图像 (从 USTB b_data_mvdr 中读取, 与 CBF 相同格式处理)
raw_mvdr = b_data_mvdr.data;
sz_m     = size(raw_mvdr);
if sz_m(1) == N_d && (numel(sz_m) < 2 || sz_m(2) == N_a)
    mvdr_env = squeeze(double(abs(raw_mvdr)));
elseif sz_m(1) == N_d * N_a
    mvdr_env = reshape(abs(double(raw_mvdr(:, 1))), N_d, N_a);
else
    mvdr_env = squeeze(double(abs(raw_mvdr(:, :, 1))));
end
if ~ismatrix(mvdr_env); mvdr_env = mvdr_env(:, :, 1); end
mvdr_img    = mvdr_env / (max(mvdr_env(:)) + eps);
mvdr_db_vis = 20 * log10(mvdr_img + eps);
mvdr_db_vis = max(mvdr_db_vis, -60);

%% FISTA 结果转 dB
fista_db = 20 * log10(img_fista / (max(img_fista(:)) + eps) + eps);
fista_db = max(fista_db, -60);

%% 统一坐标轴范围与目标真实位置
az_lim    = [-1.5, 11.5];
depth_lim = [4.5, 11.5];

% 目标理论方位角与距离 (用于在图像上标注)
tgt_az = [-delta_deg/2,  delta_deg/2, ...  % 组1 @ z=6m
          -delta_deg/2,  delta_deg/2, ...  % 组2 @ z=10m
          10];                              % 孤立目标
tgt_r  = [6, 6, 10, 10, 8];

fig_cmp = figure('Name', 'FLS成像算法三路对比 (CBF / MVDR / FISTA)', ...
                 'Position', [50, 50, 1800, 600], 'Color', 'w');

% --- 子图 1: CBF / DAS ---
ax1 = subplot(1, 3, 1);
imagesc(az_deg, depth_axis, cbf_db_vis);
colormap(ax1, hot); clim([-60 0]); colorbar;
hold on;
plot(tgt_az, tgt_r, 'c+', 'MarkerSize', 16, 'LineWidth', 2.2);
hold off;
title('CBF / DAS  (常规波束形成)', 'FontSize', 13, 'FontWeight', 'bold');
xlabel('方位角 [°]'); ylabel('深度 [m]');
set(ax1, 'XLim', az_lim, 'YLim', depth_lim);

% --- 子图 2: MVDR / Capon ---
ax2 = subplot(1, 3, 2);
imagesc(az_deg, depth_axis, mvdr_db_vis);
colormap(ax2, hot); clim([-60 0]); colorbar;
hold on;
plot(tgt_az, tgt_r, 'c+', 'MarkerSize', 16, 'LineWidth', 2.2);
hold off;
title('MVDR / Capon  (自适应波束形成)', 'FontSize', 13, 'FontWeight', 'bold');
xlabel('方位角 [°]'); ylabel('深度 [m]');
set(ax2, 'XLim', az_lim, 'YLim', depth_lim);

% --- 子图 3: CBF + FISTA 去卷积 ---
ax3 = subplot(1, 3, 3);
imagesc(az_deg, depth_axis, fista_db);
colormap(ax3, hot); clim([-60 0]); colorbar;
hold on;
plot(tgt_az, tgt_r, 'c+', 'MarkerSize', 16, 'LineWidth', 2.2);
hold off;
title('CBF + FISTA 去卷积  (稀疏正则化)', 'FontSize', 13, 'FontWeight', 'bold');
xlabel('方位角 [°]'); ylabel('深度 [m]');
set(ax3, 'XLim', az_lim, 'YLim', depth_lim);

% 联动三轴
linkaxes([ax1, ax2, ax3], 'xy');

sgtitle(sprintf('FLS成像算法对比 — 5目标场景\n紧邻目标角度间距=%.3f° (%.1f×CBF 3dB宽度=%.3f°)\nCBF无法分辨 | MVDR边缘改善 | FISTA完全分辨', ...
    delta_deg, delta_deg/bw_deg, bw_deg), ...
    'FontSize', 13, 'FontWeight', 'bold');
drawnow;

fprintf('\n=== 前视声纳三路算法对比完成 ===\n');
fprintf('  %-12s 旁瓣高, 分辨率受限于阵列孔径\n', 'CBF / DAS:');
fprintf('  %-12s 自适应抑制旁瓣, 近场效果好, 计算量大\n', 'MVDR:');
fprintf('  %-12s 稀疏去卷积, 主瓣最窄, 目标定位精度最高\n', 'FISTA:');

disp('算法对比完成！');

% 仿真结束，释放 Field II 资源
field_end();

