%% PVFC
% 

%% Inputs
% [Name]       [Description]                      [Size]
%  A            i-th body screws from i-th frame   6*n
%  M            initial relative frames M_{i,i-1}  4*4*n
%  q            joint angles                       n*1
%  qdot         joint vleocities                   n*1
%  G            link inertial matrices             6*6*n
%  Vdot_0       (optional1) base acceleration      6*1

%% Outputs
% [Name]       [Description]                      [Size]

%% Implementation
function M_q = getMassMatrix(A,M,q,G)
    %% Initialization
    n     = size(q,1);          % number of joints
     
    T    = zeros(4,4,n);        % T_{i,i-1}
    Ad_T = zeros(6,6,n);        % Ad_T_{i,i-1}

    %% A
    diagA  = zeros(6*n,n);
    for i=1:n
        diagA(6*i-5:6*i,i) = A(:,i);
    end
    
    %% G
    diagG = zeros(6*n,6*n);
    for i=1:n
        diagG(6*i-5:6*i,6*i-5:6*i) = G(:,:,i);
    end
    
    %% ad_T
    for i = 1:n
        T(:,:,i)    = exp_se3(-A(:,i)*q(i))*M(:,:,i);
        Ad_T(:,:,i) = large_Ad(T(:,:,i));
    end
    
    %% L
    L = eye(6*n);
    for i=2:n
        for j=i-1:-1:1
            L(6*i-5:6*i,6*j-5:6*j) = Ad_T(:,:,i)*L(6*i-11:6*i-6,6*j-5:6*j);
        end
    end
    
    %% Dynamic matrices
    M_q = diagA'*L'*diagG*L*diagA;
end