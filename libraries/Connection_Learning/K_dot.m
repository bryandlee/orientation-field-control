function outputArg1 = K_dot(K,dlogR_gR_inv_fro_norm,delta)
%K_DOT 이 함수의 요약 설명 위치
%   자세한 설명 위치
if nargin < 3
  delta = 1;
end
outputArg1 = K*(-1/(2*delta^2))*dlogR_gR_inv_fro_norm;
end

