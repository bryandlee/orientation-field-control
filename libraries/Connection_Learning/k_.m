function outputArg1 = k_(R_g,R,w_n,eta)
%K_ 이 함수의 요약 설명 위치
%   자세한 설명 위치
if nargin < 4
  eta = 1;
end
outputArg1 = ( 1- exp(-(norm(log_SO3(R_g'*R),'fro')^2)/(2*eta^2)) )*w_n/norm(log_SO3(R_g'*R),'fro');
end

