function outputArg1 = ddotWdepW(x,xdot,v_d,v_d_dot, Theta,index)
%DDOTWDEPW 이 함수의 요약 설명 위치
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

dPhi_11 = ddphi3D(x,1,1,index);
dPhi_12 = ddphi3D(x,1,2,index);
dPhi_13 = ddphi3D(x,1,3,index);
dPhi_22 = ddphi3D(x,2,2,index);
dPhi_23 = ddphi3D(x,2,3,index);
dPhi_33 = ddphi3D(x,3,3,index);

x_dot_1 = xdot(1);
x_dot_2 = xdot(2);
x_dot_3 = xdot(3);

W = (x_1*Theta_1 + x_2*Theta_2 + x_3*Theta_3)*Phi;

dexpW1=dexpW(W,1);
dexpW2=dexpW(W,2);
dexpW3=dexpW(W,3);

Wddot = (v_d_dot_1*Theta_1+v_d_dot_2*Theta_2+v_d_dot_3*Theta_3)*Phi ...
    +(v_d_1*Theta_1+v_d_2*Theta_2+v_d_3*Theta_3)*(x_dot_1*dPhi_1+x_dot_2*dPhi_2+x_dot_3*dPhi_3) ...
    + (x_dot_1*Theta_1+x_dot_2*Theta_2+x_dot_3*Theta_3)*(v_d_1*dPhi_1+v_d_2*dPhi_2+v_d_3*dPhi_3) ...
    + (x_1*Theta_1 + x_2*Theta_2 + x_3*Theta_3)...
    *(v_d_1*x_dot_1*dPhi_11 + v_d_1*x_dot_2*dPhi_12 + v_d_1*x_dot_3*dPhi_13 ...
    + v_d_2*x_dot_1*dPhi_12 + v_d_2*x_dot_2*dPhi_22 + v_d_2*x_dot_3*dPhi_23 ...
    + v_d_3*x_dot_1*dPhi_13 + v_d_3*x_dot_2*dPhi_23 + v_d_3*x_dot_3*dPhi_33);

Wddot_1 = Wddot(1);
Wddot_2 = Wddot(2);
Wddot_3 = Wddot(3);

outputArg1 = Wddot_1*dexpW1+Wddot_2*dexpW2+Wddot_3*dexpW3;
end

