function [t1] = GLD(x)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
k=1;
k1=0;
for k=1:1:length(x)-1
  y1=angle(x(k));
  y2=angle(x(k+1));
  if y1*y2<=0
    k1=k1+1;
  end
end
t1=k1;
end

