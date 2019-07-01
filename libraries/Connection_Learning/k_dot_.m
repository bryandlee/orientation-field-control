function outputArg1 = k_dot_(R_g,R,w_n,dlogR_gR_inv_fro_norm,eta)
%K_DOT_ 이 함수의 요약 설명 위치
%   자세한 설명 위치
if nargin < 5
  eta = 1;
end
logNorm = norm(log_SO3(R_g'*R),'fro');

outputArg1 = exp(-logNorm^2/(2*eta^2)) * 1/(2*eta^2) * dlogR_gR_inv_fro_norm * w_n/logNorm ...
    + (1-exp(-logNorm^2/(2*eta^2)))*w_n*(-1/logNorm^2) * 1/(2*logNorm) * dlogR_gR_inv_fro_norm;
end

