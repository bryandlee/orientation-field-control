function outputArg1 = w_d(R_g,R,k_,dR_g)
%W_D �� �Լ��� ��� ���� ��ġ
%   �ڼ��� ���� ��ġ
outputArg1 = k_*log_SO3(R_g*R') + dR_g*R_g';
outputArg1(1,1)=0;
outputArg1(2,2)=0;
outputArg1(3,3)=0;
end

