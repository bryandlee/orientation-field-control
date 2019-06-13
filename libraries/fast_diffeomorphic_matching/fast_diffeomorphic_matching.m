%% Fast Diffeomorphic Matching
% 2019 Bryan Dongik Lee

%% Inputs
% [Name]       [Description]                      [Size]
%  X            original data                      dim*n
%  Y            transformed data                   dim*n
%  K            num of transformations             1*1
%  mu           safety margin                      1*1
%  beta         deformation rate                   1*1

%% Outputs
% [Name]       [Description]                      [Size]
%  rho          translation variance               1*K
%  c            translation center                 dim*K
%  v            translation direction              dim*K
%  Z            transformed data                   dim*n
%  error        mean error                         1*1

%% Implementation
function [rho, c, v, Z, error] = fast_diffeomorphic_matching(X, Y, K, mu, beta)
    dim = size(X,1);

    c = zeros(dim,K);   % transformation centers
    v = zeros(dim,K);   % translation
    rho = zeros(1,K);   % variance

    options = optimoptions(@fmincon,'Algorithm', 'sqp', 'SpecifyObjectiveGradient', false, 'SpecifyConstraintGradient', false,...
        'TolCon',1e-7,'TolX',1e-10,'MaxFunEvals',1000000,'MaxIter',100000,'Display','off');

    Z = X;
    for j = 1:K
        dist = vecnorm(Z-Y);
        m = find(dist == max(dist));

        c(:,j) = Z(:,m(1));
        q = Y(:,m(1));
        v(:,j) = beta*(q - c(:,j));

        rho_max_v = exp(1/2)/(sqrt(2)*norm(v(:,j)));
        rho_0 = mu*rho_max_v/2;
        rho(:,j) = fmincon(@(rho)sum(vecnorm(locally_weighted_translation(rho,c(:,j),v(:,j),Z)-Y)),rho_0,[],[],[],[],0,mu*rho_max_v,[],options);
        Z = locally_weighted_translation(rho(:,j),c(:,j),v(:,j),Z);
    end
    error = sum(vecnorm(Z - Y))/size(X,2);
end