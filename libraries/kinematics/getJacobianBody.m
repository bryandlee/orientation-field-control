%% Body Jacobian
% 2019 Bryan Dongik Lee

%% Inputs
% [Name]      [Description]                       [Size]
%  q           joint angles                        n*1
%  A_b         screws from body frame              6*n

%% Outputs
% [Name]      [Description]                       [Size]
%  J_b         body Jacobian                       6*n

%% Examples

%% Implementation
function J_b = getJacobianBody(q,A_b)
    n = size(q,1); % number of joints
    T = eye(4,4);
    J_b = zeros(6,n);
    for i = n:-1:1
        J_b(:,i) = large_Ad(T)*A_b(:,i);
        T = T*exp_se3(-A_b(:,i)*q(i));
    end
end