%% Spatial Jacobian
% 2019 Bryan Dongik Lee

%% Inputs
% [Name]      [Description]                       [Size]
%  q           joint angles                        n*1
%  A_s         screws from body frame              6*n

%% Outputs
% [Name]      [Description]                       [Size]
%  J_s         body Jacobian                       6*n

%% Examples

%% Implementation
function J_s = getJacobianSpatial(q,A_s)
    n = size(q,1); % number of joints
    T = eye(4,4);
    J_s = zeros(6,n);
    for i = 1:n
        J_s(:,i) = large_Ad(T)*A_s(:,i);
        T = T*exp_se3(A_s(:,i)*q(i));
    end
end