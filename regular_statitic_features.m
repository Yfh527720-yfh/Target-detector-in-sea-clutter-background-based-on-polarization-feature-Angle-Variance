%对于输入的一个序列片段，计算以下时域统计特征，以 行 的形式返回特征值
function t3=regular_statitic_features(data)
% data为输入的振动信号，rpm为对应的转速信号

% 傅里叶变化获得频谱
[f,result_FFt]=transToFFT(data,20000);

% 计算频域特征

F4=RPSD(result_FFt);
t3=F4;

end

