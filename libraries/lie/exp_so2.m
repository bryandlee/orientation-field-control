function R = exp_so2(w)
%EXP_SO2 이 함수의 요약 설명 위치
%   자세한 설명 위치
    if size(w)==[2,2]
        theta = acos(w(2,1));
    elseif size(w)==[1,1]
        theta = w;
    end
    R = [cos(theta), -sin(theta) ; sin(theta),cos(theta)];
end

