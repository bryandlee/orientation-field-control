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

qi = rand(dim,1);  % init point
qf = rand(dim,1);  % final point
spline_params = rand(complexity,dim);

% constant speed reparametrization
ts = linspace(0,horizon,1e3*N);
[sp, dsp] = makeSplineP2P(qi,qf,spline_params, basis_order, horizon,ts);
spline_velocity = vecnorm(dsp);
total_length = sum(spline_velocity);

t = zeros(1,N); t(N) = horizon; j = 1;
for i = 2:N-1
    local_length = 0; k = 0;
    while true
        local_length = local_length + spline_velocity(j+k);
        if local_length > total_length/(N-1)
            t(i) = ts(j+k-1);
            j = j+k-1;
            break;
        else
            k = k+1;
        end
    end
end

Y = makeSplineP2P(qi,qf,spline_params, basis_order, horizon, t);
X = zeros(dim,N);  % flat data
for i = 1:dim
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

inv_y = Z;
for i = 1:N
    inv_y(:,i) = locally_weighted_translation_inverse(rho,c,v,inv_y(:,i));
end

%plot
[psi_X, psi_X_k] = locally_weighted_translation(rho,c,v,X);
figure(2); hold on;
axis equal; axis([0 1 0 1]);
plot(X(1,:), X(2,:),'-k','LineWidth',2);
plot(Y(1,:), Y(2,:),'--r','LineWidth',2);
plot(inv_y(1,:), inv_y(2,:),'--g','LineWidth',2);

% for i = 1:2:K
%     plot(psi_X_k(1,:,i), psi_X_k(2,:,i),'-b','LineWidth',0.1);
% end
getframe;

if(sum(vecnorm(X-inv_y))/N > 0.01)
    disp('Error: inverse mapping failed');
    return;    
end

%% Vector Field
disp('generating vector field..');

NP = 50;    % num points for vectorfield plot
x1 = linspace(0,1,NP);
x2 = x1;
x3 = [0 1];
[x1, x2, x3] = meshgrid(x1, x2, x3);


eigen_direction = (qf-qi)/norm(qf-qi);
normal_direction = null(eigen_direction(:).');
A = [eigen_direction normal_direction] * diag([-1 -5]) * [eigen_direction normal_direction]';

x1dot = A(1,1)*(x1 - qf(1)) + A(1,2)*(x2 - qf(2));
x2dot = A(2,1)*(x1 - qf(1)) + A(2,2)*(x2 - qf(2));
x3dot = - (x3 - 0);

figure(3);
subplot(1,2,1); hold on;
axis equal; axis([0 1 0 1]);
plot(X(1,:), X(2,:),'-b','LineWidth',2);
% quiver(x1, x2, x1dot, x2dot,'r');
streamslice(x1,x2,x3,x1dot,x2dot,x3dot,[],[],0);

y1 = x1;
y2 = x2;
y3 = x3;
y1dot = x1dot;
y2dot = x2dot;
y3dot = x3dot;

for l=1:NP
    disp(['mapping vector field..', num2str(100*l/NP), '%']);
    for m=1:NP
        for n=1:1

            inv_y = locally_weighted_translation_inverse(rho,c,v,[y1(l,m,n);y2(l,m,n)]);
            J = locally_weighted_translation_derivative(rho,c,v,inv_y);
            ydot = J*A*(inv_y - qf);
            
            y1dot(l,m,n) = ydot(1);
            y2dot(l,m,n) = ydot(2);
        end
    end
end

subplot(1,2,2); hold on;
axis equal; axis([0 1 0 1]);
plot(Z(1,:), Z(2,:),'-b','LineWidth',2);
% quiver(x1, x2, y1dot, y2dot ,'r');
streamslice(x1,x2,x3,y1dot,y2dot,y3dot,[],[],0);


%%
disp('done');