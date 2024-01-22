% 清除环境
clc();  % 清除命令窗口
clf();  % 清除图形窗口
clear();  % 清除工作空间变量

ncdisp("../HadISST_sst.nc")

% 设置参数
dir = 'D:/2024-MCM/2020-A/';  % 文件路径
sstname = [dir, 'HadISST_sst.nc'];  % 文件名，海表温度数据
moi = 4:7;  % 感兴趣的月份，4月到7月
latoi = [50, 75];  % 感兴趣的纬度范围（北纬）
lonoi = [-25, 10];  % 感兴趣的经度范围（东经）
landthreshold = 1/16;  % 陆地数据过滤阈值，最小NaN比例

% 读取和处理数据
sst = readnc(sstname);  % 读取海表温度数据
sst = filtbymonth(sst, moi);  % 过滤出感兴趣的月份
sst = meanbyyear(sst);  % 计算年均值
sst = trimlongitude(sst);  % 调整经度表示，使用负数表示西经
sst = filtbylatlon(sst, latoi, lonoi);  % 过滤出感兴趣的纬度和经度范围
sst.sst = permute(sst.sst, [2, 1, 3]);  % 调整数据维度顺序

% 清除陆地数据并重塑
[sst, landmask] = landclear(sst, round(landthreshold * numel(sst.time)));  % 过滤掉陆地数据
sst.sst = reshape(sst.sst, [], size(sst.sst, 3))';  % 重塑数据形状
sst.sst = sst.sst(:, any(~isnan(sst.sst)));  % 移除全为NaN的列
Y = sst.sst;  % 最终用于分析的数据


% 交叉验证和ARIMA模型选择
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
                
                % 进行预测
                [Y_pred, ~, ~] = forecast(EstMdl, length(test), 'Y0', train);
                
                % 计算RMSE
                % 这边有问题，我不知道怎么算？？？转置了维度也匹配不上
                Y_pred = Y_pred';
                rmse_values(i, j, l, k) = sqrt(mean((test - Y_pred).^2));
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


% 辅助函数定义
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

% 这里可以添加更多分析和可视化的代码
