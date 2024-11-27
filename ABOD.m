function scores = ABOD(X)  
    % ABOD: Angle-Based Outlier Detection  
    % 输入:  
    %   X - 数据矩阵，每行是一个数据点  
    %   k - 每个数据点考虑的邻居数  
    % 输出:  
    %   scores - 每个数据点的ABOD分数  
  
    n = size(X, 1); % 数据点的数量  
    scores = zeros(n, 1); % 初始化ABOD分数  
    B=1;
   
    for i = 1:n  
        X1=X-X(i,:);
        b1=1;
        k=size(X1,1);
        % 计算角度  
        angles=zeros(1,k*(k+1)/2);
        for j = 1:k  
            v1 = X(i,:) - X1(j,:);  
            v3=v1;
            v1 = v1 / norm(v1); % 单位化  
              
            for m = j+1:k  
                v2 = X(i,:) - X1(m,:); 
                v4=v2;
                v2 = v2 / norm(v2); % 单位化  
                
                % 计算两个向量之间的角度（弧度）  
                angles(b1) =acos(dot(v1, v2)/(norm(v1)*norm(v2)));  
                b1=b1+1;
            end  
        end 
        scores(B)=1/var(angles);
        B=B+1; 
        %disp(B)
        % 计算ABOD分数  
        % 这里我们简单地使用角度之和的倒数作为ABOD分数  
        % 也可以考虑其他方式，如角度之和的平方根等  
        %scores(i) = 1 / sum(angles);  
    end  
end

