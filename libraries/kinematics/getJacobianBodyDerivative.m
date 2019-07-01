%% Derivative of Body Jacobian
% 2019 Bryan Dongik Lee

%% Inputs
% [Name]      [Description]                       [Size]
%  q           joint angles                        n*1
%  A_b         screws from body frame              6*n

%% Outputs
% [Name]      [Description]                       [Size]
%  dJ_b        time derivative of body Jacobian    6*n

%% Examples

%% Implementation
function dJ = getJacobianBodyDerivative(q,A_b,qdot)
    n = size(q,1); % number of joints
    dJ = zeros(6,n);

    Ad = eye(6,6);
    for i = n:-1:1
        dJ_i = zeros(6,n);
        
        T = eye(4,4);
        
        for j = i-1:-1:1
            T = T*exp_se3(-A_b(:,j+1)*q(j+1));
            dJ_i(:,j) = large_Ad(T)*A_b(:,j);           
        end
        
        dJ = dJ + Ad * small_ad(-A_b(:,i)) * dJ_i * qdot(i);
        Ad = Ad * large_Ad(exp_se3(-A_b(:,i)*q(i)));
    end
end