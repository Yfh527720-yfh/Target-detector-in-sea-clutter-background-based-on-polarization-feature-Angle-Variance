
%频带能量
function F4=BandEnergy(f,result_FFt)
    startIndex=sum(f<=350);
    endIndex=sum(f<700);
    F4=sum(result_FFt(startIndex:endIndex));
end