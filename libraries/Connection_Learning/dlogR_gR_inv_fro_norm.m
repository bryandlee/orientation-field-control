function outputArg1= dlogR_gR_inv_fro_norm(R_g,R,R_gdot,Rdot)
%DLOGR_GR_INV_FRO_NOREM 이 함수의 요약 설명 위치
%   자세한 설명 위치
R_ginvR = R_g'*R;
R_ginvRdot = -R_g'*R_gdot*R_g'*R + R_g'*Rdot;
X = log_SO3(R_ginvR);
dX = dlog_R(R_ginvR,R_ginvRdot);
outputArg1 = 2 * trace( X' * dX );
end

