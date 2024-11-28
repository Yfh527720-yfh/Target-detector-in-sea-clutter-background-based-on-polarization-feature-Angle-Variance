clc;
clear all;
close all;
%% Read the necessary data %%

file = '19931111_163625_starea54.cdf';
ncdisp(file);
finfo = ncinfo(file);
azi = ncread(file,'azimuth_angle');
range = ncread(file,'range');
data = int16(ncread(file,'adc_data'));
N = 131072; 
n = 14;
data = permute(data,[4 3 2 1]); 
A1=8;
A2=6;
A3=11;
mode = 'auto' ;% 'raw' (no pre-processing) 
           % 'auto' (automatic pre-processing) 
           % 'dartmouth' (pre-processing for dartmouth files containing land)
pol = 'vv'; % hh, hv, vh, vv +for:VV,VH,HV,HH           
%% If try to use the data of 93, use this part to correct the data %%          
           
for a = 1:4
    for b = 1:14
        for c = 1:2
            for d = 1:131072
                if (data(d,c,b,a)) < 0
                    data(d,c,b,a) = data(d,c,b,a) + 256;
                end
            end
        end
    end
end
%% Normalization %%           
           
for rangebin = 1:n
    [I(:,rangebin),Q(:,rangebin),meanIQ,stdIQ,inbal] = NewIpixLoad(finfo,data,pol,rangebin,mode);
    R(:,rangebin) = abs(I(:,rangebin) + 1j * Q(:,rangebin));
    r(:,rangebin) = R(:,rangebin) / max(R(:,rangebin));
end
%
P1=1;
DATA=I+1j*Q;
d=32;
D=1024;
C=fix((131072-D+d)/d);
T1=zeros(1,11*C);
T2=zeros(1,11*C);
T3=zeros(1,11*C);
for b=1:1:C
  t1=GLD(DATA((d*(b-1)+1):(d*(b-1)+D),A1));
  H=hurst(R((d*(b-1)+1):(d*(b-1)+D),A1));
  t3=regular_statitic_features(R((d*(b-1)+1):(d*(b-1)+D),A1));
  if isnan(t3)==true
    t3=T3(1,P1-1);
  end
  T1(1,P1)=t1;
  T2(1,P1)=H;
  T3(1,P1)=t3;
  P1=P1+1;
end
for k=1:(A2)
  for b=1:1:C
    t1=GLD(DATA((d*(b-1)+1):(d*(b-1)+D),k));
    H=hurst(R((d*(b-1)+1):(d*(b-1)+D),k));
    t3=regular_statitic_features(R((d*(b-1)+1):(d*(b-1)+D),k));
    if isnan(t3)==true
      t3=T3(1,P1-1);
    end 
    T1(1,P1)=t1;
    T2(1,P1)=H;
    T3(1,P1)=t3;
    P1=P1+1;
  end
end
for k1=(A3):14
  for b=1:1:C
    t1=GLD(DATA((d*(b-1)+1):(d*(b-1)+D),k1));
    H=hurst(R((d*(b-1)+1):(d*(b-1)+D),k1));
    t3=regular_statitic_features(R((d*(b-1)+1):(d*(b-1)+D),k1));
    if isnan(t3)==true
      t3=T3(1,P1-1);
    end 
    T1(1,P1)=t1;
    T2(1,P1)=H;
    T3(1,P1)=t3;
    P1=P1+1;
  end
end
%
%
pol = 'vh'; % hh, hv, vh, vv +for:VV,VH,HV,HH           
%% If try to use the data of 93, use this part to correct the data %%          
           
for a = 1:4
    for b = 1:14
        for c = 1:2
            for d = 1:131072
                if (data(d,c,b,a)) < 0
                    data(d,c,b,a) = data(d,c,b,a) + 256;
                end
            end
        end
    end
end
%% Normalization %%           
           
for rangebin = 1:n
    [I(:,rangebin),Q(:,rangebin),meanIQ,stdIQ,inbal] = NewIpixLoad(finfo,data,pol,rangebin,mode);
    R(:,rangebin) = abs(I(:,rangebin) + 1i * Q(:,rangebin));
    r(:,rangebin) = R(:,rangebin) / max(R(:,rangebin));
end
%
P2=1;
DATA=I+1j*Q;
d=32;
D=1024;
C=fix((131072-D+d)/d);
T4=zeros(1,11*C);
T5=zeros(1,11*C);
T6=zeros(1,11*C);
for b=1:1:C
  t1=GLD(DATA((d*(b-1)+1):(d*(b-1)+D),A1));
  H=hurst(R((d*(b-1)+1):(d*(b-1)+D),A1));
  t3=regular_statitic_features(R((d*(b-1)+1):(d*(b-1)+D),A1));
  if isnan(t3)==true
    t3=T3(1,P1-1);
  end
  T4(1,P2)=t1;
  T5(1,P2)=H;
  T6(1,P2)=t3;
  P2=P2+1;
end
for k=1:(A2)
  for b=1:1:C
    t1=GLD(DATA((d*(b-1)+1):(d*(b-1)+D),k));
    H=hurst(R((d*(b-1)+1):(d*(b-1)+D),k));
    t3=regular_statitic_features(R((d*(b-1)+1):(d*(b-1)+D),k));
    if isnan(t3)==true
      t3=T3(1,P1-1);
    end
    T4(1,P2)=t1;
    T5(1,P2)=H;
    T6(1,P2)=t3;
    P2=P2+1;
  end
end
for k1=(A3):14
  for b=1:1:C
    t1=GLD(DATA((d*(b-1)+1):(d*(b-1)+D),k1));
    H=hurst(R((d*(b-1)+1):(d*(b-1)+D),k1));
    t3=regular_statitic_features(R((d*(b-1)+1):(d*(b-1)+D),k1));
    if isnan(t3)==true
      t3=T3(1,P1-1);
    end
    T4(1,P2)=t1;
    T5(1,P2)=H;
    T6(1,P2)=t3;
    P2=P2+1;
  end
end
%
%
pol = 'hv'; % hh, hv, vh, vv +for:VV,VH,HV,HH           
%% If try to use the data of 93, use this part to correct the data %%          
           
for a = 1:4
    for b = 1:14
        for c = 1:2
            for d = 1:131072
                if (data(d,c,b,a)) < 0
                    data(d,c,b,a) = data(d,c,b,a) + 256;
                end
            end
        end
    end
end
%% Normalization %%           
           
for rangebin = 1:n
    [I(:,rangebin),Q(:,rangebin),meanIQ,stdIQ,inbal] = NewIpixLoad(finfo,data,pol,rangebin,mode);
    R(:,rangebin) = abs(I(:,rangebin) + 1j * Q(:,rangebin));
    r(:,rangebin) = R(:,rangebin) / max(R(:,rangebin));
end
%
P3=1;
DATA=I+1j*Q;
d=32;
D=1024;
C=fix((131072-D+d)/d);
T7=zeros(1,11*C);
T8=zeros(1,11*C);
T9=zeros(1,11*C);
for b=1:1:C
  t1=GLD(DATA((d*(b-1)+1):(d*(b-1)+D),A1));
  H=hurst(R((d*(b-1)+1):(d*(b-1)+D),A1));
  t3=regular_statitic_features(R((d*(b-1)+1):(d*(b-1)+D),A1));
  if isnan(t3)==true
    t3=T3(1,P1-1);
  end 
  T7(1,P3)=t1;
  T8(1,P3)=H;
  T9(1,P3)=t3;
  P3=P3+1;
end
for k=1:(A2)
  for b=1:1:C
    t1=GLD(DATA((d*(b-1)+1):(d*(b-1)+D),k));
    H=hurst(R((d*(b-1)+1):(d*(b-1)+D),k));
    t3=regular_statitic_features(R((d*(b-1)+1):(d*(b-1)+D),k));
    if isnan(t3)==true
      t3=T3(1,P1-1);
    end
    T7(1,P3)=t1;
    T8(1,P3)=H;
    T9(1,P3)=t3;
    P3=P3+1;
  end
end
for k1=(A3):14
  for b=1:1:C
    t1=GLD(DATA((d*(b-1)+1):(d*(b-1)+D),k1));
    H=hurst(R((d*(b-1)+1):(d*(b-1)+D),k1));
    t3=regular_statitic_features(R((d*(b-1)+1):(d*(b-1)+D),k1));
    if isnan(t3)==true
      t3=T3(1,P1-1);
    end
    T7(1,P3)=t1;
    T8(1,P3)=H;
    T9(1,P3)=t3;
    P3=P3+1;
  end
end
%
%
pol = 'hh'; % hh, hv, vh, vv +for:VV,VH,HV,HH           
%% If try to use the data of 93, use this part to correct the data %%          
           
for a = 1:4
    for b = 1:14
        for c = 1:2
            for d = 1:131072
                if (data(d,c,b,a)) < 0
                    data(d,c,b,a) = data(d,c,b,a) + 256;
                end
            end
        end
    end
end
%% Normalization %%           
           
for rangebin = 1:n
    [I(:,rangebin),Q(:,rangebin),meanIQ,stdIQ,inbal] = NewIpixLoad(finfo,data,pol,rangebin,mode);
    R(:,rangebin) = abs(I(:,rangebin) + 1j * Q(:,rangebin));
    r(:,rangebin) = R(:,rangebin) / max(R(:,rangebin));
end
%
P4=1;
DATA=I+1j*Q;
d=32;
D=1024;
C=fix((131072-D+d)/d);
T10=zeros(1,11*C);
T11=zeros(1,11*C);
T12=zeros(1,11*C);
for b=1:1:C
  t1=GLD(DATA((d*(b-1)+1):(d*(b-1)+D),A1));
  H=hurst(R((d*(b-1)+1):(d*(b-1)+D),A1));
  t3=regular_statitic_features(R((d*(b-1)+1):(d*(b-1)+D),A1));
  if isnan(t3)==true
    t3=T3(1,P1-1);
  end
  T10(1,P4)=t1;
  T11(1,P4)=H;
  T12(1,P4)=t3;
  P4=P4+1;
end
for k=1:(A2)
  for b=1:1:C
    t1=GLD(DATA((d*(b-1)+1):(d*(b-1)+D),k));
    H=hurst(R((d*(b-1)+1):(d*(b-1)+D),k));
    t3=regular_statitic_features(R((d*(b-1)+1):(d*(b-1)+D),k));
    if isnan(t3)==true
      t3=T3(1,P1-1);
    end
    T10(1,P4)=t1;
    T11(1,P4)=H;
    T12(1,P4)=t3;
    P4=P4+1;
  end
end
for k1=(A3):14
  for b=1:1:C
    t1=GLD(DATA((d*(b-1)+1):(d*(b-1)+D),k1));
    H=hurst(R((d*(b-1)+1):(d*(b-1)+D),k1));
    t3=regular_statitic_features(R((d*(b-1)+1):(d*(b-1)+D),k1));
    if isnan(t3)==true
      t3=T3(1,P1-1);
    end
    T10(1,P4)=t1;
    T11(1,P4)=H;
    T12(1,P4)=t3;
    P4=P4+1;
  end
end
j2=0;
%
T=zeros(11*C,12);
T=[T1;T2;T3;T4;T5;T6;T7;T8;T9;T10;T11;T12]';
T13=T(2*C+1:9*C,1:12);
T14=T(1:C,1:12);
T15=T(C+1:2*C,1:12);
k = 5; % 选择邻居数  
scores = ABOD(T13); % 计算ABOD分数  

% 排序分数并找到异常点  
[U, idx] = sort(scores);  
%outliers = U(end-20:end); % 假设我们关注最后10个最高分的数据点作为异常点 
%disp(idx(end-50:end-21))
cankao=T13(idx(end-(10+fix(8*C/100)):end-(fix(8*C/100))));%虚警率可控
disp(cankao)
for j=1:8*C
   ang(j)=ABODb(T13(j,:) ,cankao);
end
%disp(ang)
[so, idx1] = sort(ang);
z=so(end-(12+fix(8*C/100)));
%disp(z)
parfor j1=1:C
   ang1(j1)=ABODb(T14(j1,:) ,cankao);
   %disp(ang1(j1))
   if ang1(j1)>z
       j2=j2+1;
   end
end
j4=0;
parfor j3=1:1:C
   ang1(j3)=ABODb(T15(j3,:) ,cankao);
   %disp(ang1(j1))
   if ang1(j3)<=z
       j2=j2+1;
   end
end
%disp(ang1)
disp((j2+j4)/(2*C))
