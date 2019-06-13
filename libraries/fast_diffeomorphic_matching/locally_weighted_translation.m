%% Locally Weighted Translation
% 2019 Bryan Dongik Lee

%% Inputs
% [Name]       [Description]                      [Size]
%  rho          translation variance               1*K
%  c            translation center                 dim*K
%  v            translation direction              dim*K
%  X            data                               dim*n

%% Outputs
% [Name]       [Description]                      [Size]
%  psi_X        transformed data                   dim*n
%  psi_X_k      (optional)all k-transformed data   dim*n*K

%% Implementation
function [psi_X, psi_X_k] = locally_weighted_translation(rho, c, v, X)
    K = size(rho,2);
    [dim, n] = size(X);
    psi_X_k = zeros(dim,n,K);
    
    psi_X = X;
    for k = 1:K
        psi_X = psi_X + v(:,k)*exp(-rho(k)^2 * vecnorm(psi_X-c(:,k)).^2);
        psi_X_k(:,:,k) = psi_X;
    end
end