function R = exp_so2(w)
%EXP_SO2 �� �Լ��� ��� ���� ��ġ
%   �ڼ��� ���� ��ġ
    if size(w)==[2,2]
        theta = acos(w(2,1));
    elseif size(w)==[1,1]
        theta = w;
    end
    R = [cos(theta), -sin(theta) ; sin(theta),cos(theta)];
end

