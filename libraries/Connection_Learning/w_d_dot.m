function outputArg1 = w_d_dot(x,xdot,R,Rdot,Theta,index, v_d,v_d_dot, k_, k_dot_,R_g)
%W_D_DOT 이 함수의 요약 설명 위치
%   자세한 설명 위치
R_gRInv = R_g*R';
logR_gRInv = log_SO3(R_gRInv);
R_g_dot = dR_g(x,xdot,Theta,index);
R_gRInv_dot = R_g_dot*R' - R_g*R'*Rdot*R';
dtdlogR_gRInv = dlog_R(R_gRInv,R_gRInv_dot);

First = k_dot_*logR_gRInv;
Second = k_*dtdlogR_gRInv;

[outputArg1,outputArg2] = dtdv_ddxexpW(x,xdot,v_d,v_d_dot, Theta,index);

Third =  outputArg1 * R_g';
Fourth = outputArg2 * R_g';
Fifth = - dR_g(x,v_d,Theta,index)  * R_g' * R_g_dot * R_g';

outputArg1 = First+ Second + Third + Fourth + Fifth;

outputArg1(1,1)= 0;
outputArg1(2,2)= 0;
outputArg1(3,3)= 0;

end

