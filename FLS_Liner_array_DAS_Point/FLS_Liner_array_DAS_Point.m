%% 完整工作流程
%{
1. 场景构建
   ├─ 读取 STL 模型 (ball.stl)
   ├─ 体素化下采样
   ├─ 旋转、缩放、平移
   └─ 生成海底散射体

2. 信号设计
   ├─ Chirp 脉冲
   ├─ 匹配滤波器
   └─ 窗函数加权

3. Field II 仿真
   ├─ 创建阵列 (128元)
   ├─ 设置发射聚焦
   ├─ 计算散射响应
   └─ 匹配滤波 + 加窗

4. 波束成形
   ├─ DAS 算法
   ├─ 扇形扫描
   └─ 孔径加权

5. 成像显示
   └─ 70 dB 动态范围图像
%}
%% 添加路径
clear all; close all;
addpath(genpath("./ustb"))
addpath(genpath("./Field_II"))
% addpath(genpath("./Lib"))
%% 基本常量
c0 = 1500;  % 声速 (m/s)
fs = 10e6;  % 采样频率 10 MHz
dt = 1/fs;  % 采样间隔 0.1 微秒

field_init(0);
set_field('c',c0);
set_field('fs',fs);
set_field('use_rectangles',1);

% 探头参数，线阵，x轴放置，z为前视方向
probe = uff.linear_array();
f0                      = 0.9e+06;      % 中心频率 900 kHz
lambda                  = c0/f0;           % 波长 λ = c/f = 1.67 mm
probe.element_height    = 15e-3;           % 阵元高度 15 mm
probe.pitch             = lambda;        % 阵元间距 1.67 mm
kerf                    = 0.1*lambda;        % 阵元间隙 0.167 mm
probe.element_width     = probe.pitch-kerf;     % 实际宽度
% lens_el                 = 20e-3;           % 仰角焦点的位置
probe.N                 = 128;            % 128个阵元
pulse_duration          = 2.5;             % 脉冲持续时间 [周期]

%% 用户参数 
stl_file = 'D:\109_Research\Matlab_code\FLS_Liner_array_DAS_Point\ball.stl';  % STL 模型
stl_scale   = 0.05;                     % 🔶 比例因子 0.5=缩小一半
stl_pos     = [2 8 6];                % 📌 放置位置 (x ,z ,y[depth])，注意 y 被用于 depth
stl_rot_deg = [0 0 90];               % 欧拉角，[滚转, 俯仰, 偏航]
%海底地形生成
seabed_len  = 15;  seabed_wid = 15;    % 📌 z×x 尺寸   m
seabed_y0   = 10;  seabed_amp = 3;     % 📌 平均深度 / 起伏（depth 在 y 方向）
N_seabed_scatter = 1e4;                % 海底散射点数

FOV_yaw   = -45:1:45;                  % 🔶 水平角网格 91 条
FOV_pitch = -45:1:45;                  % 🔶 俯仰角网格 91 条
FOV_R     = seabed_len;                % 射线长度

%% ① STL 处理
try
    TR = stlread(stl_file); stl_V = TR.Points;
catch
    [~, stl_V] = stlread(stl_file);
end
% 去重
[stl_V, ~, idxFaces] = unique(stl_V, 'rows');   % ~ 去重后点数可能更小

% 体素化下采样
voxel = 0.4;               % 体素边长（与 STL 单位同量级）
minCoord = floor(min(stl_V)/voxel);   %找到每列点的最小坐标,除以体素大小并向下取整
idx = round((stl_V - minCoord*voxel)/voxel);   % 映射到体素索引

[~, firstIdx] = unique(idx, 'rows', 'stable');
stl_V = stl_V(firstIdx, :);

% 旋转和缩放
stl_V = stl_scale * stl_V;             % 🔶 缩放
r = stl_rot_deg;
Rxm = [1 0 0;0 cosd(r(1)) -sind(r(1));0 sind(r(1)) cosd(r(1))]; %绕x轴旋转
Rym = [cosd(r(2)) 0 sind(r(2));0 1 0;-sind(r(2)) 0 cosd(r(2))]; %绕y轴旋转
Rzm = [cosd(r(3)) -sind(r(3)) 0;sind(r(3)) cosd(r(3)) 0;0 0 1]; %绕z轴旋转
Rot  = Rzm*Rym*Rxm; %组合旋转，注意顺序

% 📌 坐标变换并重排为 (x, y[depth], z)
pos_target = (Rot * stl_V.').' + stl_pos;  % 原为 Nx3 的 [x z y]       .' 表示转置
pos_target = pos_target(:, [1 3 2]);       % 改为 [x y z]，其中 y 为 depth   坐标重排  ：表示取所有，这里是所有行


%% ② 海底散射体生成（以 y 为深度）
% 生成网格
Nx = round( sqrt(N_seabed_scatter * seabed_wid / seabed_len) ); %round 四舍五入函数
Nz = round( N_seabed_scatter / Nx );

x_lin = linspace(-seabed_wid/2, seabed_wid/2, Nx);   % 横向 x 100个点，linspace 等差序列生成函数
z_lin = linspace(0, seabed_len, Nz);                % 纵向 z 100个点（“前方”）
[Xg,Zg] = meshgrid(x_lin, z_lin);

%加入起伏
Yg = seabed_y0 + seabed_amp * ( ...
      0.5*sin(2*pi*Xg/seabed_wid) + 0.5*sin(2*pi*Zg/seabed_len) ) + ...
      0.1*randn(size(Xg));  % 📌 y 方向作为深度，加入扰动

pos_seabed = [Xg(:) Yg(:) Zg(:)];   %网格展开为点列表

% 设置散射系数
amp_target = 0.23*ones(size(pos_target,1),1);           %目标均匀反射
amp_seabed = 0.2 * abs(randn(size(pos_seabed,1),1));    %海底不均匀散射
%% ③ 场景预览 (3-D + 两视图 + 网格)
figure('Name','FLS Scene Preview','Position',[50 60 1400 700]);

% --- 3-D ---------------------------------------------------------------
subplot(2,3,1);
scatter3(pos_seabed(:,1),pos_seabed(:,2),pos_seabed(:,3),3,'b.'); hold on;
scatter3(pos_target(:,1),pos_target(:,2),pos_target(:,3),10,'g','filled');
plot3(0,0,0,'rp','MarkerFaceColor','r','MarkerSize',10);
set(gca,'ZDir','reverse'); axis equal tight; grid on;
xlabel('x');ylabel('y');zlabel('z'); title('3-D');
view(35,25);

% --- 顶视 (x-y) --------------------------------------------------------
subplot(2,3,2);
scatter(pos_seabed(:,1),pos_seabed(:,2),3,'b.'); hold on;
scatter(pos_target(:,1),pos_target(:,2),10,'g','filled');
plot(0,0,'rp','MarkerFaceColor','r','MarkerSize',8);   
axis equal tight; grid on; xlabel('x');ylabel('y'); title('Top x-y');

% --- 侧视 (y-z) --------------------------------------------------------
subplot(2,3,3);
scatter(pos_seabed(:,2),pos_seabed(:,3),3,'b.'); hold on;
scatter(pos_target(:,2),pos_target(:,3),10,'g','filled');
plot(0,0,'rp','MarkerFaceColor','r','MarkerSize',8);   
set(gca,'YDir','reverse'); axis equal tight; grid on;
xlabel('y');ylabel('z'); title('Side y-z');

% --- FOV 网格 ----------------------------------------------------------
subplot(2,3,[4 5 6]); hold on; axis equal; grid on;
set(gca,'ZDir','reverse'); xlabel('x');ylabel('y');zlabel('z');
title('FOV Grid & Scatterers');

% 水平网格 (pitch=0)
for ya = FOV_yaw
    dir = [ sind(ya)  cosd(ya)  0 ] * FOV_R;
    plot3([0 dir(1)],[0 0],[0 dir(2)],'b-');
end
% 垂直网格 (yaw=0)
for pa = FOV_pitch
    dir = [ 0  cosd(pa)  sind(pa)] * FOV_R;
    plot3([0 dir(1)],[0 dir(3)],[0 dir(2)],'r-');
end
% 绘制稀疏海底和目标
scatter3(pos_seabed(1:15:end,1),pos_seabed(1:15:end,2), ...
         pos_seabed(1:15:end,3),3,'b.');
scatter3(pos_target(:,1),pos_target(:,2),pos_target(:,3),10,'g','filled');
plot3(0,0,0,'rp','MarkerFaceColor','r','MarkerSize',10);
legend('Yaw grid','Pitch grid','scatterers','target','array','Location','best');
view(3);
drawnow;


%% 激励脉冲生成
pulse = uff.pulse();                    %创建一个uff标准的脉冲
pulse.fractional_bandwidth = 1/9;       % 11%带宽
pulse.center_frequency = f0;            % 中心频率 900kHz
BW=pulse.fractional_bandwidth*f0;       % 带宽 100kHz
pulse_len=2e-3;                         % 2 毫秒 = 0.002 秒

% Tx = xdc_piston(0.25*lambda,0.25*lambda/4);

% Chirp信号 线性调频
% t0 = (-1/pulse.fractional_bandwidth/f0): dt : (1/pulse.fractional_bandwidth/f0);    %时间轴定义
impulse_response = 1;       % 理想冲激响应
t_p=0:1/fs:pulse_len;       % 生成时间向量
excitation = chirp(t_p,f0-BW/2,pulse_len,f0+BW/2).*hamming(numel(t_p))';  % 加Hamming窗减少旁瓣  % 逐元素乘法  .*
excitation = excitation-mean(excitation); % 去除直流分量
one_way_ir = conv(impulse_response,excitation); % 卷积生成单程脉冲

%匹配滤波使信噪比最大化
MF = conj(flipud(excitation(:)));                  % FIR 形式, 时间反转 + 共轭
two_way_ir = conv(one_way_ir,MF);                  % 脉冲压缩
lag = length(two_way_ir)/2+1;                       % 调整峰值位置
% 我们显示脉冲以检查延迟（lag）估计是否准确
% ）并且检查脉冲是否对称）
figure;
plot((0:(length(two_way_ir)-1))*dt -lag*dt,two_way_ir); hold on; grid on; axis tight
plot((0:(length(two_way_ir)-1))*dt -lag*dt,abs(hilbert(two_way_ir)),'r')
plot([0 0],[min(two_way_ir) max(two_way_ir)],'g');
legend('2-ways pulse','Envelope','Estimated lag');
title('2-ways impulse response Field II');

% 构建aperture，阵列创建
% 子阵元剖分，更精确的声场计算
noSubAz=round(probe.element_width/(lambda/8));        % 子阵元数 （方位角方向）
noSubEl=round(probe.element_height/(lambda/8));       % 子阵元数 （仰角方向）
Th = xdc_linear_array (probe.N, probe.element_width, probe.element_height, kerf, noSubAz, noSubEl, [0 0 Inf]);  % 发射阵列

% Th = xdc_piston(lambda/80, lambda/800);
Rh = xdc_linear_array (probe.N, probe.element_width, probe.element_height, kerf, noSubAz, noSubEl, [0 0 Inf]);  % 接收阵列，与发射相同

% 我们还将激励、脉冲响应和挡板设置如下
xdc_excitation (Th, excitation);        % 设置之前生成的 chirp 激励信号
xdc_impulse (Th, impulse_response);     % 设置换能器的冲激响应
xdc_baffle(Th, 0);                      % 设置阵列的障板条件 - `0` = 软障板（soft baffle）
xdc_center_focus(Th,[0 0 0]);           % 设置聚焦中心位置  - `[0 0 0]` = 坐标原点（阵列中心）
xdc_impulse (Rh, impulse_response);
xdc_baffle(Rh, 0);
xdc_center_focus(Rh,[0 0 0]);

%% 声场计算
% 平面波序列
F_number = 1;                               % F-number = 焦距 / 孔径宽度
alpha_max = 0;                              %不进行角度扫描
Na=1;                                      % 平面波数量，只发射一个平面波 
F=1;                                        % 帧数
alpha=linspace(-alpha_max,alpha_max,Na);    % 角度向量

% 定义散射体模型
point_position=[pos_target;pos_seabed]; 
point_amplitudes=[amp_target;amp_seabed];

% 输出数据
cropat=round(2*50e-3/c0/dt);    % 计算时间窗长度
CPW=zeros(cropat,probe.N,1,1);  % 初始化数据矩阵，存储原始的通道数据 CPW = Compounded Plane Wave（复合平面波）

% Compute CPW signals
time_index=0;
disp('Field II: Computing CPW dataset');
wb = waitbar(0, 'Field II: Computing Forward Looking MBES dataset');
win = blackman(probe.N); %布莱克曼窗
for f=1:F
    waitbar(0, wb, sprintf('Field II: Computing CPW dataset, frame %d',f));
    for n=1:Na
        waitbar(n/Na, wb);
        disp(['Calculating angle ',num2str(n),' of ',num2str(Na)]);
         
        % 设置发射孔径
        tmp = zeros(1,probe.N/2);
        xdc_apodization(Th,0,[tmp(1:end-1),1,1,tmp(1:end-1)]);
        xdc_times_focus(Th,0,probe.geometry(:,1)'.*sin(alpha(n))/c0);  %发射聚焦延迟
        
        % 设置接收孔径
        xdc_apodization(Rh, 0, ones(1,probe.N));
        xdc_focus_times(Rh, 0, zeros(1,probe.N)); %接收不聚焦

        % 计算散射响应
        [v,t]=calc_scat_multi(Th, Rh, point_position, point_amplitudes);

        % 匹配滤波+加窗
        % 降低旁瓣，改善图片对比度
        sig_mf = zeros(size(v));
        for m = 1:probe.N
            tmp = conv(v(:,m), MF, 'same');
            sig_mf(:,m) = tmp*win(m);
        end
        % 存储数据
        CPW(1:size(v,1),:,n,f)=sig_mf;
        
        %
        clear sig_mf tmp
        % Save transmit sequence
        seq(n)=uff.wave();
        seq(n).probe=probe;
        % seq(n).source.azimuth=alpha(n);
        seq(n).source.distance=Inf;
        seq(n).sound_speed=c0;
        seq(n).delay = -lag*dt+t;
    end
end
close(wb);
%% 数据存储（UFF格式）
% Channel Data
% 
% 在代码的这一部分，我们创建一个 UFF 数据结构，专门用于存储捕获的超声波通道数据。
channel_data = uff.channel_data();
channel_data.sampling_frequency = fs;
channel_data.sound_speed = c0;
channel_data.initial_time = 0;
channel_data.pulse = pulse;
channel_data.probe = probe;
channel_data.sequence = seq;
channel_data.data = CPW./max(CPW(:)); % 归一化

%% 扫描定义
depth_axis=linspace(0e0,25e0,512).';        % 深度：0-25m，512点
azimuth_axis=linspace(-45,45,451)/180*pi;   % 方位：±45°，451点


sca=uff.sector_scan('azimuth_axis',azimuth_axis,'depth_axis',depth_axis);
%% 波束成形（DAS算法）
%
% 有了 channel_data（通道数据）和 scan（扫描定义），我们就拥有了生成超声/声纳图像所需的全部信息。
% 现在我们使用 USTB 的 pipeline（处理流水线）结构，它除了需要channel_data 和 scan 之外,
% 还需要一个 apodization（孔径加权）结构。

pipe=pipeline();
pipe.channel_data=channel_data;
pipe.scan=sca;

pipe.receive_apodization.window=uff.window.tukey25;
pipe.receive_apodization.f_number=F_number;

%% 成像显示
%
% *pipeline* 结构允许你通过组合多个内置的 *processes* 来实现不同的波束成形器。
% 通过改变 *process* 链，可以实现其他波束成形序列。 
% 它会返回另一个*UFF* 结构：*beamformed_data*（波束成形后的数据）。
% 
% 为了实现本例的目标，我们使用延迟求和算法（在 *das_mex()* 进程中实现）以及相干复合。 

b_data=pipe.go({midprocess.das()});     % postprocess.coherent_compounding()});

% Display images
b_data.plot();
colormap(hot); 
colorbar; 
caxis([-70 0]); %动态范围 70 dB

