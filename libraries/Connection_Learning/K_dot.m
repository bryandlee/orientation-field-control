function outputArg1 = K_dot(K,dlogR_gR_inv_fro_norm,delta)
%K_DOT �� �Լ��� ��� ���� ��ġ
%   �ڼ��� ���� ��ġ
if nargin < 3
  delta = 1;
end
outputArg1 = K*(-1/(2*delta^2))*dlogR_gR_inv_fro_norm;
end

