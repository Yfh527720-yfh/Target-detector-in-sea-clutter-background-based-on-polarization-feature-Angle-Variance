
function [f,result_FFt]=transToFFT(data,fs)
%用于求解数据的快速傅里叶变换结果，以便于绘制频谱图
N=length(data);
%去均值
data=data-mean(data);
%求频率分辨率
df=fs/(N-1);
f=(0:N-1)*df;
Y=fft(data)/N*2;
result_FFt=abs(Y);
%取一半
result_FFt=result_FFt(1:ceil(N/2));
f=f(1:ceil(N/2));
end
