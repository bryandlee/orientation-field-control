%% Locally Weighted Translation
% 2019 Bryan Dongik Lee

%% Inputs
% [Name]       [Description]                      [Size]
%  rho          translation variance               1*1
%  c            translation center                 dim*1
%  v            translation direction              dim*1
%  Z            transformed data                   dim*1

%% Outputs
% [Name]       [Description]                      [Size]
%  X            original data                      dim*1

%% Implementation
function X = locally_weighted_translation_inverse(rho, c, v, Z)
    count = 0;
    max_iter = 1e5;
    max_error = 1e-10;
    stepsize = 0.2;
    
    X = Z;

    K = size(rho,2);
    for k = K:-1:1
        while true
            [psi_X, psi_X_k] = locally_weighted_translation(rho(k),c(:,k),v(:,k),X);

            % Newton's method
    %         X = X - (locally_weighted_translation_derivative(rho,c,v,X,psi_X_k))\(psi_X - Z);

            % gradient descent
            X = X - stepsize*2*locally_weighted_translation_derivative(rho(k),c(:,k),v(:,k),X,psi_X_k)'*(psi_X - Z);

            if norm(locally_weighted_translation(rho(k),c(:,k),v(:,k),X) - Z) < max_error
                Z = X;
                break;
            end

            count = count + 1;
            if count > max_iter
                disp('locally_weighted_translation_inverse: max count exceeded');
                return;
            end
        end
    end 
end