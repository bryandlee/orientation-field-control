function omega = log_SO3( R )
%LOG_SO3 이 함수의 요약 설명 위치
%   자세한 설명 위치
s = size(R);
dim = s(1); 
R=reshape(R,[dim,dim]);
if dim == 2
    a = acos(R(1,1));
    b = -a;
    if sin(a)-(R(2,1)) < 1.0e-10
        theta = a;
    elseif sin(b)-(R(2,1)) < 1.0e-10
        theta = b;
    end
    omega = [0,-theta;theta,0];
elseif dim == 3
    if trace(R)==-1
        stop
    elseif trace(R)==3
        omega = zeros(3);
    else
        theta = acos((trace(R)-1)/2);
        omega = theta/(2*sin(theta)) * (R-R');
    end
end

end

