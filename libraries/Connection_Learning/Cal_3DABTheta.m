
function [A,B,Theta] = Cal_3DABTheta(T,W,index)
%CAL_ABTHETA 이 함수의 요약 설명 위치
%   자세한 설명 위치
s = size(T);
N = s(3);
dim = 3;
x = T(1:dim,dim+1,:);
R = T(1:dim,1:dim,:);
ss = size(index);
M = ss(1);

A = zeros(dim,M,dim*(dim-1)/2);
B = zeros(dim,dim,M,M);
r = zeros(N,dim*(dim-1)/2);

for i=1:N
    r(i,:) = skew(reshape(log_SO3(R(:,:,i)),[dim,dim]));
end

for i=1:N
    for j=1:dim
        A(j,:,:) = A(j,:,:) + reshape(x(j,1,i) * phi3D(reshape(x(:,1,i),[dim,1]),index) * r(i,:),[1,M,dim*(dim-1)/2]);
    end
    for j=1:dim
        for k=1:dim
            B(j,k,:,:) = B(j,k,:,:) + reshape( x(j,1,i) * x(k,1,i)* phi3D(reshape(x(:,1,i),[dim,1]),index) * phi3D(reshape(x(:,1,i),[dim,1]),index)', [1,1,M,M] ) ;
        end
    end
end

if dim==2
    Theta_together = ([W+reshape(B(1,1,:,:),[M,M]), reshape(B(1,2,:,:),[M,M]); 
        reshape(B(2,1,:,:),[M,M]), W+reshape(B(2,2,:,:),[M,M])])\([ reshape(A(1,:,:),[M,1]) ; reshape(A(2,:,:),[M,1]) ]);
    Theta = zeros(dim,1,M);
    Theta(1,:,:) = Theta_together(1:M,:)';
    Theta(2,:,:) = Theta_together(M+1:2*M,:)';
elseif dim==3
    Theta_together = ([W+reshape(B(1,1,:,:),[M,M]), reshape(B(1,2,:,:),[M,M]), reshape(B(1,3,:,:),[M,M]); 
        reshape(B(2,1,:,:),[M,M]), W+reshape(B(2,2,:,:),[M,M]), reshape(B(2,3,:,:),[M,M]); 
        reshape(B(3,1,:,:),[M,M]), reshape(B(3,2,:,:),[M,M]), W+reshape(B(3,3,:,:),[M,M])])\([reshape(A(1,:,:),[M,3]) ; reshape(A(2,:,:),[M,3]) ; reshape(A(3,:,:),[M,3])]);
    Theta = zeros(dim,dim,M);
    Theta(1,:,:) = Theta_together(1:M,:)';
    Theta(2,:,:) = Theta_together(M+1:2*M,:)';
    Theta(3,:,:) = Theta_together(2*M+1:3*M,:)';
else
    pass
end

end

