function H = hurst(time_series)
    % 步骤1: 计算时序数据的均值标准差
    n = length(time_series);
    mean_ts = mean(time_series);
    std_ts = std(time_series);
    % 步骤2: 构建新的时间序列数组
    new_ts = cumsum(time_series - mean_ts);
    
    % 步骤3: 计算不同尺度的范围标准差
    range = zeros(ceil(log2(n)), 1);
    scale = zeros(ceil(log2(n)), 1);
    
    count = 1;
    for i = 16:16:n
        range(count) = max(new_ts(i-15:i)) - min(new_ts(i-15:i));
        scale(count) = log2(i);
        count = count + 1;
    end
    
    % 步骤4: 对范围标准差进行线性拟合
    p = polyfit(scale,log(range),1);
    
    % 步骤5: 计算Hurst指数
    H = p(1);
end
