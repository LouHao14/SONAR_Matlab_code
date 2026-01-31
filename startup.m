%% Matlab 项目环境初始化脚本
% 作用：自动识别当前项目根目录，并将其子文件夹加入 Matlab 路径

% 1. 获取当前脚本所在的绝对路径（即项目根目录）
rootPath = fileparts(mfilename('fullpath'));

% 2. 添加 common_utils 及其所有子文件夹
% genpath 会递归生成该目录下所有子目录的路径列表
addpath(genpath(fullfile(rootPath, 'common_utils')));

% 3. 添加 FLS 和 MBES 的源码路径（可选，建议添加）
addpath(genpath(fullfile(rootPath, 'FLS', 'src')));
addpath(genpath(fullfile(rootPath, 'MBES', 'src')));

% 4. 打印提示信息
fprintf('正在初始化声呐仿真项目...\n');
fprintf('项目根目录: %s\n', rootPath);
fprintf('已成功加载 common_utils (ustb, Field_II) 及相关组件。\n');