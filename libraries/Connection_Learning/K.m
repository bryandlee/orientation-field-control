function outputArg1 = K(R_g,R,delta)
%K �� �Լ��� ��� ���� ��ġ
%   �ڼ��� ���� ��ġ
if nargin < 3
  delta = 1;
end
outputArg1 = exp(-(norm(log_SO3(R_g'*R),'fro')^2)/(2*delta^2));
end

