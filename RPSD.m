
%相对功率谱熵
function F5=RPSD(result_FFt)

    Fi=result_FFt.*result_FFt;
    Pi=Fi/sum(Fi);
    F5=-sum(Pi.*log2(Pi))/log2(length(result_FFt));
end
