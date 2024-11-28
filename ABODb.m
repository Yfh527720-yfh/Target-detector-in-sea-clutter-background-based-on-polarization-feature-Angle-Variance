function score = ABODb(X,Y)  
    % ABOD: Angle-Based Outlier Detection  
    % 输入:  
    %   X - 数据矩阵，每行是一个数据点  
    %   k - 每个数据点考虑的邻居数  
    % 输出:  
    %   scores - 每个数据点的ABOD分数  
  
        [k,n]= size(Y); % 数据点的数量  
        score = zeros(k, 1); % 初始化ABOD分数  
        %disp(k)
  
        % 计算角度 
        b=k*(k-1)/2; 
        b1=1;
        B=1;
        angles = zeros(b, 1);  
        for j = 1:k  
            v1 = X - Y(j,:);  
            v3=v1; 
            v1 = v1 / norm(v1); % 单位化  
           
            for m = j+1:k  
                v2 = X - Y(m,:);  
                v4=v2;
                v2 = v2 / norm(v2); % 单位化  
                  
                % 计算两个向量之间的角度（弧度）  
                angles(b1) =acos(dot(v1, v2)/(norm(v1)*norm(v2)));  
                b1=b1+1;  
            end  
            
        end  
        score=1/var(angles);
        % 计算ABOD分数  
        % 这里我们简单地使用角度之和的倒数作为ABOD分数  
        % 也可以考虑其他方式，如角度之和的平方根等  
        %scores(i) = 1 / sum(angles);  
    
end

