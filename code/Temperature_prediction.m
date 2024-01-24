%% 清除环境
clc();  % 清除命令窗口
clear();  % 清除工作空间变量

ncdisp("../HadISST_sst.nc")

%% 设置参数
dir = 'D:/2024-MCM/2020-A/';  % 文件路径
sstname = [dir, 'HadISST_sst.nc'];  % 文件名，海表温度数据
moi = 4:7;  % 感兴趣的月份，4月到7月
latoi = [50, 70];  % 感兴趣的纬度范围（北纬）
lonoi = [-30, 10];  % 感兴趣的经度范围（东经）
landthreshold = 1/16;  % 陆地数据过滤阈值，最小NaN比例

%% 读取和处理数据
sst = readnc(sstname);  % 读取海表温度数据
sst = filtbymonth(sst, moi);  % 过滤出感兴趣的月份
sst = meanbyyear(sst);  % 计算年均值
sst = trimlongitude(sst);  % 调整经度表示，使用负数表示西经
sst = filtbylatlon(sst, latoi, lonoi);  % 过滤出感兴趣的纬度和经度范围
sst.sst = permute(sst.sst, [2, 1, 3]);  % 调整数据维度顺序

% 基本统计信息
disp('基本统计信息:');
disp(['最小值: ', num2str(min(sst.sst(:)))]);

% 数据形状
disp('数据形状:');
disp(size(sst.sst));

%% 绘图
figure;
% 设置坐标轴的实际刻度
x = linspace(-30, 10, size(sst.sst, 2));  % 假设数据的经度范围是 -30 到 10
y = linspace(50, 70, size(sst.sst, 1));   % 假设数据的纬度范围是 50 到 70
% 绘制数据
data = flipud(sst.sst(:, :, 154));
imagesc(x, y, data);
set(gca, 'YDir', 'normal');  % 确保Y轴是从下到上的
colorbar;
title('SST 分布（第一个时间点）');
xlabel('经度');
ylabel('纬度');

% 自定义颜色映射
colormap([0.5 0.5 0.5; jet]);  % 在jet颜色映射前添加灰色，对应于NaN值

% 保留原有颜色映射的同时，将NaN值设置为灰色
nanMask = isnan(sst.sst(:, :, 1));

% 绘制海岸线
hold on;
load coastlines;  % 加载海岸线数据
plot(coastlon, coastlat, 'k-', 'LineWidth', 0.7);  % 使用加载的海岸线数据绘制海岸线
hold off;

% 设置坐标轴范围和刻度
axis ([-30, 10, 50, 70]);  % 调整坐标轴范围
set(gca, 'XTick', -30:10:10);  % 设置x轴刻度
set(gca, 'YTick', 50:5:70);   % 设置y轴刻度
grid on;  % 显示网格







%% 检查特定数据点
% 例如，查看特定位置和时间的温度
disp('特定数据点:');
disp(sst.sst(18, 30, 1)); % 示例：第10纬度，第10经度，第一个时间点的温度

% 初始化一个数组来存储每个时间点的缺失值数量
missingValues = zeros(1, 154);

% 遍历时间点1到154
for timePoint = 1:154
    % 计算每个时间点的缺失值数量
    missingValues(timePoint) = sum(isnan(sst.sst(:, :, timePoint)), 'all');
end

% 初始化一个标志变量来跟踪是否存在变化
hasChange = false;

% 检查并输出缺失值数量的变化
for timePoint = 2:154
    if missingValues(timePoint) ~= missingValues(timePoint - 1)
        disp(['时间点 ', num2str(timePoint), ' 的缺失值数量从 ', ...
              num2str(missingValues(timePoint - 1)), ' 变化到 ', ...
              num2str(missingValues(timePoint))]);
        hasChange = true;
    end
end

% 如果没有变化，则输出相应信息
if ~hasChange
    disp('从时间点 1 到 154 的缺失值数量始终无变化。');
end

% 初始化两个数组，一个用于索引和时间点，另一个用于温度
indexTimeChanges = []; % 格式为 [latIndex, lonIndex, timePoint]
temperatureChanges = []; % 格式为 [beforeTemp, afterTemp]，均为字符串

% 遍历时间点，从第二个时间点开始
for timePoint = 2:size(sst.sst, 3)
    % 遍历每个经纬度位置
    for latIndex = 1:size(sst.sst, 1)
        for lonIndex = 1:size(sst.sst, 2)
            % 获取当前和前一个时间点的温度
            currentTemp = sst.sst(latIndex, lonIndex, timePoint);
            previousTemp = sst.sst(latIndex, lonIndex, timePoint - 1);

            % 检查海水变海冰或海冰变海水的情况
            if (currentTemp < -273.15 && previousTemp >= -273.15) || (currentTemp >= -273.15 && previousTemp < -273.15)
                % 记录索引和时间点
                indexTimeChanges = [indexTimeChanges; latIndex, lonIndex, timePoint];
                % 格式化并记录温度变化
                temperatureChanges = [temperatureChanges; {sprintf('%.2f', previousTemp), sprintf('%.2f', currentTemp)}];
            end
        end
    end
end

% 检查是否有变化并输出
if isempty(indexTimeChanges)
    disp('没有发现海水与海冰之间的状态变化。');
else
    disp('发现海水与海冰之间的状态变化：');
    for i = 1:size(indexTimeChanges, 1)
        fprintf('经度索引: %d, 纬度索引: %d, 时间点: %d, 前温度: %s°C, 后温度: %s°C\n', ...
                indexTimeChanges(i, 1), indexTimeChanges(i, 2), indexTimeChanges(i, 3), ...
                temperatureChanges{i, 1}, temperatureChanges{i, 2});
    end
end

%% 清除陆地数据并重塑
% [sst, landmask] = landclear(sst, round(landthreshold * numel(sst.time)));  % 过滤掉陆地数据
% sst.sst = reshape(sst.sst, [], size(sst.sst, 3))';  % 重塑数据形状
% sst.sst = sst.sst(:, any(~isnan(sst.sst)));  % 移除全为NaN的列
Y = sst.sst;  % 最终用于分析的数据


%% 交叉验证和ARIMA模型选择
% 定义交叉验证的折数
K = 6;

% 随机生成交叉验证的索引
n = size(Y, 1);
indices = crossvalIndices(n, K);

% 定义要测试的ARIMA模型参数
p_values = [0, 1, 2];  % 自回归项
d_values = [0, 1, 2];  % 差分项
q_values = [0, 1, 2];  % 移动平均项

% 初始化存储RMSE的矩阵
rmse_values = zeros(length(p_values), length(d_values), length(q_values), K);

for k = 1:K
    % 分割数据为训练集和测试集
    train = Y(indices ~= k, :);
    test = Y(indices == k, :);

    for i = 1:length(p_values)
        for j = 1:length(d_values)
            for l = 1:length(q_values)
               % 拟合ARIMA模型
                mdl = arima(p_values(i), d_values(j), q_values(l));
                train = train(:);
                EstMdl = estimate(mdl, train);
                
                % 获取预测值，限制预测长度
                Y_pred = forecast(EstMdl, length(test), 'Y0', train);
                length_test = length(test);
                size_Y_pred = size(Y_pred, 1);
                disp(['Length of test: ', num2str(length_test)]);
                disp(['Size of Y_pred (dimension 1): ', num2str(size_Y_pred)]);

                % 对齐维度
                disp(['Size of test: ', num2str(size(test))]);
                disp(['Size of Y_pred: ', num2str(size(Y_pred))]);
                aligned_test = test(1:min(size(test, 1), numel(Y_pred)), :);
                aligned_test = aligned_test(:);  % 将 aligned_test 转换为列向量
                aligned_test = aligned_test(1:length(Y_pred));
                disp(['Size of aligned_test: ', num2str(size(aligned_test))]);


                % 输出关键变量尺寸
                disp(['Size of aligned_test: ', num2str(size(aligned_test))]);
                disp(['Size of Y_pred: ', num2str(size(Y_pred))]);

                % 计算RMSE
                rmse_values(i, j, l, k) = sqrt(mean((aligned_test - Y_pred).^2, 'all'));
            end
        end
    end
end

% 计算平均RMSE
mean_rmse = mean(rmse_values, 4);

% 选择RMSE最小的模型参数
[min_rmse, idx] = min(mean_rmse(:));
[p_idx, d_idx, q_idx] = ind2sub(size(mean_rmse), idx);
best_p = p_values(p_idx);
best_d = d_values(d_idx);
best_q = q_values(q_idx);

% 输出最佳模型参数
fprintf('Best ARIMA Model: p=%d, d=%d, q=%d\n', best_p, best_d, best_q);


%% 辅助函数定义
function indices = crossvalIndices(totalSize, numGroups)
    % 生成交叉验证的随机索引
    indices = randperm(totalSize);
    indices = mod(indices, numGroups) + 1;
end

% 读取NetCDF文件的函数
function data = readnc(filename)
    data.sst = ncread(filename, 'sst');  % 读取海表温度数据，大小为 360x180x1847
        
    time_data_bnds = ncread(filename, 'time_bnds'); % 读取'time'变量数据
    start = datetime(1870, 1, 1);
    time = start + days(time_data_bnds(1, :));
    data.time = time';
    
    data.lat = ncread(filename, 'latitude');  % 读取纬度数据
    data.lon = ncread(filename, 'longitude');  % 读取经度数据
end

function data = filtbymonth(data, months)
    dates = data.time;
    monthVec = month(dates);
    monthIdx = ismember(monthVec, months);                        
    data.sst = data.sst(:, :, monthIdx);
    data.time = data.time(monthIdx);
end

function data = meanbyyear(data)
    dates = data.time;
    years = year(dates);
    uniqueYears = unique(years);
    meanSST = NaN(size(data.sst, 1), size(data.sst, 2), length(uniqueYears));
    for i = 1:length(uniqueYears)
        yearIdx = years == uniqueYears(i);
        meanSST(:, :, i) = mean(data.sst(:, :, yearIdx), 3, 'omitnan');
    end
    data.sst = meanSST;
    data.time = uniqueYears;
end

function data = trimlongitude(data)
    % 将经度从[0, 360]调整为[-180, 180]
    shiftedLon = data.lon > 180;
    data.lon(shiftedLon) = data.lon(shiftedLon) - 360;

    % 对经度进行排序，并获取排序后的索引
    [data.lon, sortIdx] = sort(data.lon);

    % 重新排序sst数据
    % 注意: 根据您提供的信息，sst的维度顺序是 [longitude, latitude, time]
    data.sst = data.sst(sortIdx, :, :);
end


% 根据经纬度过滤数据的函数
function data = filtbylatlon(data, latrange, lonrange)
    % 找到符合范围的经纬度索引
    latIdx = data.lat >= latrange(1) & data.lat <= latrange(2);
    lonIdx = data.lon >= lonrange(1) & data.lon <= lonrange(2);

    % 过滤数据
    % 注意: 根据您的数据，sst的维度顺序是 [longitude, latitude, time]
    data.sst = data.sst(lonIdx, latIdx, :);
    data.lat = data.lat(latIdx);
    data.lon = data.lon(lonIdx);
end


function [data, landmask] = landclear(data, threshold)
    landmask = sum(isnan(data.sst), 3) / size(data.sst, 3) > threshold;
    data.sst(repmat(landmask, [1, 1, size(data.sst, 3)])) = NaN;
end

%% 这里可以添加更多分析和可视化的代码
