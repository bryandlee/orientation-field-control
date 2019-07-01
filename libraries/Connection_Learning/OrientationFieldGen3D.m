function R = OrientationFieldGen3D(x,Theta,index)
%UNTITLED 이 함수의 요약 설명 위치
%   자세한 설명 위치
ss = size(index);
N = ss(1);
x_1 = x(1);
x_2 = x(2);
x_3 = x(3);
Theta_1 = reshape(Theta(1,:,:),[3,N]);
Theta_2 = reshape(Theta(2,:,:),[3,N]);
Theta_3 = reshape(Theta(3,:,:),[3,N]);
Phi = phi3D(x,index);

a = (x_1*Theta_1 + x_2*Theta_2 + x_3*Theta_3)*Phi;

R = exp_so3(skew(a));
end

