function r = FKJXGX(x, y)  
    % 检查输入是否为复数向量  
    
  
    % 获取序列长度  
    N = length(x);  
    M = length(y);  
      
    % 填充零以匹配较长的序列长度  
    if N < M  
        x = [zeros(1, M-N) x];  
    else  
        y = [zeros(1, N-M) y];  
    end  
      
    % 初始化互相关结果  
    r = zeros(1, N + M - 1);  
      
    % 计算互相关系数  
    for lag = -(N + M - 2) : (N + M - 2)  
        if lag >= 0  
            x_lag = x(lag + 1 : end);  
        else  
            x_lag = [zeros(1, -lag) x];  
        end  
          
        r(lag + N) = sum(conj(x_lag) .* y);  
    end  
      
    % 归一化（可选）  
    r = r / (norm(x) * norm(y));  
end  
  


