function outputArg1 = v_ddtdR_g(x,xdot,v_d,Theta,index)
%V_DDTDR_G 이 함수의 요약 설명 위치
%   자세한 설명 위치
ss = size(index);
N = ss(1);
Theta_1 = reshape(Theta(1,:,:),[3,N]);
Theta_2 = reshape(Theta(2,:,:),[3,N]);
Theta_3 = reshape(Theta(3,:,:),[3,N]);

x_1 = x(1);
x_2 = x(2);
x_3 = x(3);

v_d_1 = v_d(1);
v_d_2 = v_d(2);
v_d_3 = v_d(3);

Phi = phi3D(x,index);

dPhi_1 = dphi3D(x,1,index);
dPhi_2 = dphi3D(x,2,index);
dPhi_3 = dphi3D(x,3,index);

x_dot_1 = xdot(1);
x_dot_2 = xdot(2);
x_dot_3 = xdot(3);

W = (x_1*Theta_1 + x_2*Theta_2 + x_3*Theta_3)*Phi;

W_dot_v_d = (v_d_1*Theta_1 + v_d_2*Theta_2 + v_d_3*Theta_3)*Phi ...
    + (x_1*Theta_1 + x_2*Theta_2 + x_3*Theta_3)*(v_d_1*dPhi_1 + v_d_2*dPhi_2 + v_d_3*dPhi_3);
W_dot_x_dot = (x_dot_1*Theta_1 + x_dot_2*Theta_2 + x_dot_3*Theta_3)*Phi ...
    + (x_1*Theta_1 + x_2*Theta_2 + x_3*Theta_3)*(x_dot_1*dPhi_1 + x_dot_2*dPhi_2 + x_dot_3*dPhi_3);

outputArg1 = ddexpW(W,1,1) * W_dot_v_d(1) * W_dot_x_dot (1) ...
    + ddexpW(W,1,2) * W_dot_v_d(1) * W_dot_x_dot (2) ...
    + ddexpW(W,1,3) * W_dot_v_d(1) * W_dot_x_dot (3) ...
    + ddexpW(W,1,2) * W_dot_v_d(2) * W_dot_x_dot (1) ...
    + ddexpW(W,2,2) * W_dot_v_d(2) * W_dot_x_dot (2) ...
    + ddexpW(W,2,3) * W_dot_v_d(2) * W_dot_x_dot (3) ...
    + ddexpW(W,1,3) * W_dot_v_d(3) * W_dot_x_dot (1) ...
    + ddexpW(W,2,3) * W_dot_v_d(3) * W_dot_x_dot (2) ...
    + ddexpW(W,3,3) * W_dot_v_d(3) * W_dot_x_dot (3);

end

