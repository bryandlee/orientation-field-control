%% Derivative of Locally Weighted Translation
% 2019 Bryan Dongik Lee

%% Inputs
% [Name]       [Description]                      [Size]
%  rho          translation variance               1*K
%  c            translation center                 dim*K
%  v            translation direction              dim*K
%  X            data                               dim*1
%  psi_X_k      (optional)all k-transformed data   dim*1*K

%% Outputs
% [Name]       [Description]                      [Size]
%  dpsi_X       derivative at X                    dim*dim

%% Implementation
function dpsi_X = locally_weighted_translation_derivative(rho, c, v, X, varargin)
    if nargin < 5
        [~, psi_X_k] = locally_weighted_translation(rho,c,v,X);
    else
        psi_X_k = varargin{1};
    end
    
    K = size(rho,2);
    dim = size(X,1);
    
    dpsi_X = eye(dim) - 2*(rho(1)^2)*exp(-rho(1)^2 * norm(X-c(:,1))^2) * v(:,1) * (X-c(:,1))';
    for k=2:K
        dpsi_X = (eye(dim) - 2*(rho(k)^2)*exp(-rho(k)^2 * norm(psi_X_k(:,k-1)-c(:,k))^2)*v(:,k)*(psi_X_k(:,k-1)-c(:,k))') * dpsi_X;
    end
end