function outputArg1 = K(R_g,R,delta)
%K 이 함수의 요약 설명 위치
%   자세한 설명 위치
if nargin < 3
  delta = 1;
end
outputArg1 = exp(-(norm(log_SO3(R_g'*R),'fro')^2)/(2*delta^2));
end

