function outputArg1 = dR_g(x,v,Theta,index)
%DR_G 이 함수의 요약 설명 위치
%   자세한 설명 위치
ss = size(index);
N = ss(1);

Theta_1 = reshape(Theta(1,:,:),[3,N]);
Theta_2 = reshape(Theta(2,:,:),[3,N]);
Theta_3 = reshape(Theta(3,:,:),[3,N]);

Phi = phi3D(x,index);

dPhi_1 = dphi3D(x,1,index);
dPhi_2 = dphi3D(x,2,index);
dPhi_3 = dphi3D(x,3,index);

v_1 = v(1);
v_2 = v(2);
v_3 = v(3);
x_1 = x(1);
x_2 = x(2);
x_3 = x(3);

W = (x_1*Theta_1+x_2*Theta_2+x_3*Theta_3)*Phi;

First_term = (v_1*Theta_1 + v_2*Theta_2 + v_3*Theta_3)*Phi;
Second_term = (x_1*Theta_1 + x_2*Theta_2 + x_3*Theta_3)*(v_1*dPhi_1 + v_2*dPhi_2 + v_3*dPhi_3);

Wdot = First_term+Second_term;

Wdot_1 = Wdot(1);
Wdot_2 = Wdot(2);
Wdot_3 = Wdot(3);

dexpW1=dexpW(W,1);
dexpW2=dexpW(W,2);
dexpW3=dexpW(W,3);

outputArg1 = dexpW1*Wdot_1 + dexpW2*Wdot_2 + dexpW3*Wdot_3;
end

