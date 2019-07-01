function outputArg1 = dexpW(w,axis)
%DEXPW 이 함수의 요약 설명 위치
%   자세한 설명 위치
if size(w) == [3,3]
    w = skew(w);
end
normw = norm(w);
w_1 = w(1);
w_2 = w(2);
w_3 = w(3);
dw_1 = ([0,0,0;0,0,-1;0,1,0]);
dw_2 = ([0,0,1;0,0,0;-1,0,0]);
dw_3 = ([0,-1,0;1,0,0;0,0,0]);
cw = cos(normw);
sw = sin(normw);
skew_w = skew(w);
if normw < 1.0e-02
    if axis == 1
       outputArg1 = (-1/2)*w_1*skew_w + dw_1 ...
           + (-1/12)*w_1*skew_w^2 + (1/2) * (dw_1*skew_w+skew_w*dw_1);
    elseif axis ==2
       outputArg1 = (-1/2)*w_2*skew_w + dw_2 ...
           + (-1/12)*w_2*skew_w^2 + (1/2) * (dw_2*skew_w+skew_w*dw_2); 
    elseif axis ==3
       outputArg1 = (-1/2)*w_3*skew_w + dw_3 ...
           + (-1/12)*w_3*skew_w^2 + (1/2) * (dw_3*skew_w+skew_w*dw_3); 
    end
else
    if axis == 1
       outputArg1 = ((cw/normw^2)-(sw/normw^3))*w_1*skew_w + (sw/normw)*dw_1 ...
           + ((sw/normw^3)-2*((1-cw)/normw^4))*w_1*skew_w^2 + ((1-cw)/(normw^2)) * (dw_1*skew_w+skew_w*dw_1);
    elseif axis ==2
       outputArg1 = ((cw/normw^2)-(sw/normw^3))*w_2*skew_w + (sw/normw)*dw_2 ...
           + ((sw/normw^3)-2*((1-cw)/normw^4))*w_2*skew_w^2 + ((1-cw)/(normw^2)) * (dw_2*skew_w+skew_w*dw_2); 
    elseif axis ==3
       outputArg1 = ((cw/normw^2)-(sw/normw^3))*w_3*skew_w + (sw/normw)*dw_3 ...
           + ((sw/normw^3)-2*((1-cw)/normw^4))*w_3*skew_w^2 + ((1-cw)/(normw^2)) * (dw_3*skew_w+skew_w*dw_3); 
    end
end
end

