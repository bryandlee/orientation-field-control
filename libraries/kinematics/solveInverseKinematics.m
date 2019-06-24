%% Forward Kinematics Solver
% 2018 Bryan Dongik Lee

%% Inputs
% [Name]     [Description]                       [Size]
%  A_s        screws from spatial frame           6*n
%  T_s        desired transformation              4*4
%  q0         (optional) inital q                 n*1

%% Outputs
% [Name]     [Description]                       [Size]
%  q          joint angles                        n*1

%% Examples

%% Implementation
function q = solveInverseKinematics(A_s, T_s, varargin)
    EPS = 1e-10;
    maxtries = 1e4;
    
    n = size(A_s,2); % number of joints
    
    if nargin > 2
        q = varargin{1};
        if size(q) ~= [n,1]
            error('Inverse Kinematics Error: wrong input size');
        end
    else
        q = zeros(n,1);
    end
    
    tries = 0;
    while true
        T = eye(4,4);
        for i = 1:n
            T = T*exp_se3(A_s(:,i)*q(i));
        end

        S = log_SE3(T_s * inverse_SE3(T));
        norm(S)
        if norm(S) < EPS
            break
        end
        
        J = getJacobianSpatial(q, A_s);
        q = q + 0.01*pinv(J)*S;
        
        tries = tries+1;
        if tries > maxtries
            disp('Inverse Kinematics Error: Max Try Exceeded');
            break
        end
    end
    
    q = wrapToPi(q);
end