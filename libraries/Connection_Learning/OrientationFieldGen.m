function R = OrientationFieldGen( x,Theta,order_phi )
%UNTITLED 이 함수의 요약 설명 위치
%   자세한 설명 위치
s = size(x);
dim = s(1);
a = zeros(dim*(dim-1)/2,1);
ss = size(Phi(x,order_phi));
M = ss(1);
for i=1:dim
    a = a + reshape(x(i)*reshape(Theta(i,:,:),[dim*(dim-1)/2,M])*Phi(x,order_phi),[dim*(dim-1)/2,1]);
end
if dim==3
    R = exp_so3(skew(a));
elseif dim==2
    R = exp_so2((a));
end

