clc(); clf(); clear(); 
import('preprocess.*'); 

% 参数设置
dir = 'D:/2024-MCM/2020-A/'; 
sstname = [dir,'HadISST_sst.nc']; 

moi = 4:7; % 感兴趣的月份
latoi = [50,75]; % 感兴趣的纬度范围（北纬）
lonoi = [-25,10]; % 感兴趣的经度范围（东经）

geodetic2cartesian(0,0,latoi,lonoi,1); % 将地理坐标转换为笛卡尔坐标
landthreshold = 1/16; % 陆地区域的NaN值比例阈值
sst = readnc(sstname); % 读取.nc格式的海表温度数据文件
sst.time = datetime(1800,1,1) + sst.time; % 将时间转换为从1800年1月1日开始的日期格式
sst = filtbymonth(sst, moi); % 按月份过滤数据
sst = meanbyyear(sst); % 计算每年的平均值
sst = trimlongitude(sst); % 调整经度表示，使用西经表示而非180度以上
sst = filtbylatlon(sst, latoi, lonoi); % 按纬度和经度过滤数据
sst.sst = permute(sst.sst, [2, 1, 3]); % 调整数据维度，方便后续处理
[sst, landmask] = landclear(sst, round(landthreshold * numel(sst.time))); % 过滤掉陆地区域
sst.sst = reshape(sst.sst, [], 60).'; % 重塑数据形状
sst.sst = sst.sst(:, any(sst.sst)); % 去除全为NaN的列
Y = sst.sst; % 最终的数据集

% 数据分析和模型选择
A = [
    crossvalidation(Y, [1:50, []], 51:60)
    crossvalidation(Y, [1:40, 51:60], 41:50)
    crossvalidation(Y, [1:30, 41:60], 31:40)
    crossvalidation(Y, [1:20, 31:60], 21:30)
    crossvalidation(Y, [1:10, 21:60], 11:20)
    crossvalidation(Y, [ [], 11:60], 1:10)
];
A = sort(sqrt(A)); % 计算并排序RMSE（均方根误差）

% 绘图
b = boxplot(A, 'Labels', {'(p=1, d=1)', '(p=1, d=2)', '(p=2, d=1)', '(p=2, d=2)'});
grid('on'); % 开启网格
set(gca, 'YScale', 'Log'); % 设置Y轴为对数刻度
title(['Choose ARIMA model by using RMSE', newline, 'with six fold cross validation']); % 标题
xlabel('Parameters'); % X轴标签
ylabel('Root Mean Square Error( ˛a)'); % Y轴标签
set(gca, 'fontsize', 22, 'fontname', 'times new roman'); % 设置字体大小和类型

% 下面是辅助函数的定义，包括交叉验证、模型拟合和预测等
