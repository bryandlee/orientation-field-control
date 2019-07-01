function outputArg1 = w_d(R_g,R,k_,dR_g)
%W_D 이 함수의 요약 설명 위치
%   자세한 설명 위치
outputArg1 = k_*log_SO3(R_g*R') + dR_g*R_g';
outputArg1(1,1)=0;
outputArg1(2,2)=0;
outputArg1(3,3)=0;
end

