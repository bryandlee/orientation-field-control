function outputArg1 = dlog_R(R,Rdot)
%DLOG_R 이 함수의 요약 설명 위치
%   자세한 설명 위치
if trace(R)==-1
    stop
elseif trace(R)==3
    logR = zeros(3);
    theta = 0;
else
    theta = acos((trace(R)-1)/2);
    logR = theta/(2*sin(theta)) * (R-R');
end
if theta < 1.0e-8
    outputArg1 = -1/8*trace(Rdot)*(R-R') + 1/2*(Rdot+R'*Rdot*R');
else
    outputArg1 = trace(Rdot)*( theta*(cos(theta))/(sin(theta)) -1 )*(R-R')/(4*sin(theta)^2) + (theta)/(2*sin(theta))*(Rdot+R'*Rdot*R');
end

