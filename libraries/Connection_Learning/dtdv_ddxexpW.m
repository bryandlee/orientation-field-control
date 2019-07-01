function [outputArg1,outputArg2] = dtdv_ddxexpW(x,xdot,v_d,v_d_dot, Theta,index)
%DTDV_DDXEXPW 이 함수의 요약 설명 위치
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

v_d_dot_1 = v_d_dot(1);
v_d_dot_2 = v_d_dot(2);
v_d_dot_3 = v_d_dot(3);

Phi = phi3D(x,index);

dPhi_1 = dphi3D(x,1,index);
dPhi_2 = dphi3D(x,2,index);
dPhi_3 = dphi3D(x,3,index);

x_dot_1 = xdot(1);
x_dot_2 = xdot(2);
x_dot_3 = xdot(3);

W = (x_1*Theta_1 + x_2*Theta_2 + x_3*Theta_3)*Phi;

dexpW1=dexpW(W,1);
dexpW2=dexpW(W,2);
dexpW3=dexpW(W,3);

ddexpW11 =  ddexpW(W,1,1);
ddexpW12 =  ddexpW(W,1,2);
ddexpW13 =  ddexpW(W,1,3);
ddexpW22 =  ddexpW(W,2,2);
ddexpW23 =  ddexpW(W,2,3);
ddexpW33 =  ddexpW(W,3,3);

dPhi_11 = ddphi3D(x,1,1,index);
dPhi_12 = ddphi3D(x,1,2,index);
dPhi_13 = ddphi3D(x,1,3,index);
dPhi_22 = ddphi3D(x,2,2,index);
dPhi_23 = ddphi3D(x,2,3,index);
dPhi_33 = ddphi3D(x,3,3,index);

W_dot = (v_d_dot_1*Theta_1+v_d_dot_2*Theta_2+v_d_dot_3*Theta_3)*Phi ...
    + (v_d_1*Theta_1+v_d_2*Theta_2+v_d_3*Theta_3) * (x_dot_1*dPhi_1+x_dot_2*dPhi_2+x_dot_3*dPhi_3) ...
    + (x_dot_1*Theta_1+x_dot_2*Theta_2+x_dot_3*Theta_3) * (v_d_1*dPhi_1+v_d_2*dPhi_2+v_d_3*dPhi_3) ...
    + (x_1*Theta_1 + x_2*Theta_2 + x_3*Theta_3) * (v_d_dot_1*dPhi_1+v_d_dot_2*dPhi_2+v_d_dot_3*dPhi_3) ...
    + (x_1*Theta_1 + x_2*Theta_2 + x_3*Theta_3) ... 
    * (v_d_1 * x_dot_1 * dPhi_11 + v_d_1 * x_dot_2 * dPhi_12 + v_d_1 * x_dot_3 * dPhi_13 ...
    + v_d_2 * x_dot_1 * dPhi_12 + v_d_2 * x_dot_2 * dPhi_22 + v_d_2 * x_dot_3 * dPhi_23 ...
    + v_d_3 * x_dot_1 * dPhi_13 + v_d_3 * x_dot_2 * dPhi_23 + v_d_3 * x_dot_3 * dPhi_33); 

outputArg1 = W_dot(1)*dexpW1 + W_dot(2)*dexpW2 + W_dot(3)*dexpW3;

W_ddot1 = (v_d_1*Theta_1+v_d_2*Theta_2+v_d_3*Theta_3)*Phi ...
    + (x_1*Theta_1 + x_2*Theta_2 + x_3*Theta_3) * (v_d_1*dPhi_1+v_d_2*dPhi_2+v_d_3*dPhi_3);
W_ddot2 = (x_dot_1*Theta_1+x_dot_2*Theta_2+x_dot_3*Theta_3)*Phi ...
    + (x_1*Theta_1 + x_2*Theta_2 + x_3*Theta_3) * (x_dot_1*dPhi_1+x_dot_2*dPhi_2+x_dot_3*dPhi_3);

outputArg2 = W_ddot1(1) * W_ddot2(1) * ddexpW11 + W_ddot1(1) * W_ddot2(2) * ddexpW12 + W_ddot1(1) * W_ddot2(3) * ddexpW13 ...
    + W_ddot1(2) * W_ddot2(1) * ddexpW12 + W_ddot1(2) * W_ddot2(2) * ddexpW22 + W_ddot1(2) * W_ddot2(3) * ddexpW23 ...
    + W_ddot1(3) * W_ddot2(1) * ddexpW13 + W_ddot1(3) * W_ddot2(2) * ddexpW23 + W_ddot1(3) * W_ddot2(3) * ddexpW33;
end

