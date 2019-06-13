close all
clear 
clc

%% Settings
dim = 2;    % data dimension
N = 100;    % num of data

K = 150;    % num of transformations
mu = 0.9;   % safety margin
beta = 0.5; % deformation rate

%% Data generation
disp('generating data..');

% data settings
complexity  = 3;   % num of spline coefs per joints
basis_order = 4;   % cubic spline
horizon  = 1;      % trajectory horizon

Y = zeros(dim,N);  % curve data
X = zeros(dim,N);  % flat data
t = linspace(0,horizon,N);

qi = rand(dim,1);  % init point
qf = rand(dim,1);  % final point
spline_params = rand(complexity,dim);
p2ptrajectory = makeSplineP2P(qi,qf,spline_params, basis_order, horizon);

for i = 1:dim
    Y(i,:) = fnval(p2ptrajectory(i), t);
    X(i,:) = linspace(Y(i,1),Y(i,N),N);
end


% plot
figure(1);
subplot(1,2,1); hold on;
axis equal; axis([0 1 0 1]);
plot(X(1,:), X(2,:),'-b','LineWidth',2);
plot(Y(1,:), Y(2,:),'--k','LineWidth',2);

%% Diffeomorphic Matching
disp('diffeomorphic matching..');

[rho, c, v, Z, matching_error] = fast_diffeomorphic_matching(X,Y,K,mu,beta);
if matching_error > 0.01
    disp('Error: diffeomorphic matching failed. try different trajectory');
    return;
end

% plot
subplot(1,2,2); hold on;
axis equal; axis([0 1 0 1]);
plot(Z(1,:), Z(2,:),'-b','LineWidth',2);
plot(Y(1,:), Y(2,:),'--k','LineWidth',2);
getframe;

%% inverse map
disp('verifying inverse map..');

Z_inv = Z;
for i = 1:N
    for k=K:-1:1
        Z_inv(:,i) = locally_weighted_translation_inverse(rho(k),c(:,k),v(:,k),Z_inv(:,i));
    end   
end

%plot
[psi_X, psi_X_k] = locally_weighted_translation(rho,c,v,X);
figure(2); hold on;
axis equal; axis([0 1 0 1]);
plot(X(1,:), X(2,:),'-k','LineWidth',2);
plot(Y(1,:), Y(2,:),'--r','LineWidth',2);
plot(Z_inv(1,:), Z_inv(2,:),'--g','LineWidth',2);

% for i = 1:2:K
%     plot(psi_X_k(1,:,i), psi_X_k(2,:,i),'-b','LineWidth',0.1);
% end
getframe;

if(sum(vecnorm(X-Z_inv))/N > 0.01)
    disp('Error: inverse mapping failed');
    return;    
end

%% Vector Field
disp('generating vector field..');

NP = 20;    % num points for vectorfield plot
x1 = linspace(0,1,NP);
x2 = x1;
x3 = [0 1];
[x1, x2, x3] = meshgrid(x1, x2, x3);


eigen_direction = (qf-qi)/norm(qf-qi);
normal_direction = [eigen_direction(2); -eigen_direction(1)];
A = [eigen_direction normal_direction]' * [-1/5 0; 0 -4/5] * [eigen_direction normal_direction];

x1dot = A(1,1)*(x1 - qf(1)) + A(1,2)*(x2 - qf(2));
x2dot = A(2,1)*(x1 - qf(1)) + A(2,2)*(x2 - qf(2));
x3dot = - (x3 - 0);

figure(3);
subplot(1,2,1); hold on;
axis equal; axis([0 1 0 1]);
plot(X(1,:), X(2,:),'-b','LineWidth',2);
quiver(x1, x2, x1dot, x2dot,'r');
streamslice(x1,x2,x3,x1dot,x2dot,x3dot,[],[],0);

x1dot_trans = x1dot;
x2dot_trans = x2dot;
x3dot_trans = x3dot;

for l=1:NP
    disp(['mapping vector field..', num2str(100*l/NP), '%']);
    for m=1:NP
        for n=1:2
            Z_inv = [x1(l,m,n);x2(l,m,n)];
            for k=K:-1:1
                Z_inv = locally_weighted_translation_inverse(rho(k),c(:,k),v(:,k),Z_inv);
            end
            J = locally_weighted_translation_derivative(rho,c,v,Z_inv);
            xdot_trans = J*[x1dot(l,m,n);x2dot(l,m,n)];
            xdot_trans_norm = norm(xdot_trans);
            
            x1dot_trans(l,m,n) = xdot_trans(1)/xdot_trans_norm;
            x2dot_trans(l,m,n) = xdot_trans(2)/xdot_trans_norm;
        end
    end
end

subplot(1,2,2); hold on;
axis equal; axis([0 1 0 1]);
plot(Z(1,:), Z(2,:),'-b','LineWidth',2);
quiver(x1, x2, x1dot_trans, x2dot_trans ,'r');
streamslice(x1,x2,x3,x1dot_trans,x2dot_trans,x3dot,[],[],0);


%%
disp('done');