close all
clear 
clc

%% Settings
dim = 3;    % data dimension
N = 100;    % num of data

K = 150;    % num of transformations
mu = 0.9;   % safety margin
beta = 0.5; % deformation rate

%% Data Generation
robot = makeFrankaPanda;
qi = getFeasiblePose(robot);
qf = getFeasiblePose(robot);

for i = 1:robot.dof
    q(i,:) = linspace(qi(i), qf(i), 100);
end
for i = 1:100
    T_d(:,:,i) = solveForwardKinematics(q(:,i),robot.A, robot.M);
end

visualizePanda(q);

%% Data Load
% load('../data/3.mat');
Y = reshape(T_d(1:3,4,:),[dim N]);
X = zeros(dim,N);  % flat data
for i = 1:dim
    X(i,:) = linspace(Y(i,1),Y(i,N),N);
end

qi = Y(:,1);
qf = Y(:,end);

plot_min = min(Y,[],2) - 0.2;
plot_max = max(Y,[],2) + 0.2;

% plot
figure(2);
subplot(1,2,1); hold on;
axis equal; axis([plot_min(1) plot_max(1) plot_min(2) plot_max(2) plot_min(3) plot_max(3)]);
view([135 35]); xlabel('x'); ylabel('y'); zlabel('z');

plot3(X(1,:), X(2,:), X(3,:),'-b','LineWidth',2);
plot3(Y(1,:), Y(2,:), Y(3,:), '--k','LineWidth',2);

%% Diffeomorphic Matching
disp('diffeomorphic matching..');

[rho, c, v, Z, matching_error] = fast_diffeomorphic_matching(X,Y,K,mu,beta);
if matching_error > 0.01
    disp('Error: diffeomorphic matching failed. try different trajectory');
    return;
end

% plot
subplot(1,2,2); hold on;
axis equal; axis([plot_min(1) plot_max(1) plot_min(2) plot_max(2) plot_min(3) plot_max(3)]);
view([135 35]); xlabel('x'); ylabel('y'); zlabel('z');

plot3(Z(1,:), Z(2,:), Z(3,:),'-b','LineWidth',2);
plot3(Y(1,:), Y(2,:), Y(3,:), '--k','LineWidth',2);
getframe;

%% inverse map
disp('verifying inverse map..');

inv_y = Z;
for i = 1:N
    inv_y(:,i) = locally_weighted_translation_inverse(rho,c,v,inv_y(:,i));
end

%plot
[psi_X, psi_X_k] = locally_weighted_translation(rho,c,v,X);
figure(3); hold on;
axis equal; axis([plot_min(1) plot_max(1) plot_min(2) plot_max(2) plot_min(3) plot_max(3)]);
view([135 35]); xlabel('x'); ylabel('y'); zlabel('z');

plot3(X(1,:), X(2,:), X(3,:), '-k','LineWidth',2);
plot3(Y(1,:), Y(2,:), Y(3,:),'--r','LineWidth',2);
plot3(inv_y(1,:), inv_y(2,:), inv_y(3,:),'--g','LineWidth',2);

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

NP = 2;    % num points for vectorfield plot
x1 = linspace(plot_min(1),plot_max(1),NP);
x2 = linspace(plot_min(2),plot_max(2),NP);
x3 = linspace(plot_min(3),plot_max(3),NP);
[x1, x2, x3] = meshgrid(x1, x2, x3);

eigen_direction = (qf-qi)/norm(qf-qi);
normal_direction = null(eigen_direction(:).');
A = [eigen_direction normal_direction] * diag([-1 -5 -5]) * [eigen_direction normal_direction]';

x1dot = A(1,1)*(x1 - qf(1)) + A(1,2)*(x2 - qf(2)) + A(1,3)*(x3 - qf(3));
x2dot = A(2,1)*(x1 - qf(1)) + A(2,2)*(x2 - qf(2)) + A(2,3)*(x3 - qf(3));
x3dot = A(3,1)*(x1 - qf(1)) + A(3,2)*(x2 - qf(2)) + A(3,3)*(x3 - qf(3));

figure(4);
subplot(1,2,1); hold on;
axis equal; axis([plot_min(1) plot_max(1) plot_min(2) plot_max(2) plot_min(3) plot_max(3)]);
view([135 35]); xlabel('x'); ylabel('y'); zlabel('z');

plot3(X(1,:), X(2,:), X(3,:), '-r','LineWidth',2);

for i = 1:100
    startpoints(:,i) = (plot_max - plot_min) .* rand(3,1) + plot_min;
end
streamline(x1,x2,x3,x1dot,x2dot,x3dot,startpoints(1,:),startpoints(2,:),startpoints(3,:));

y1 = x1;
y2 = x2;
y3 = x3;
y1dot = x1dot;
y2dot = x2dot;
y3dot = x3dot;

for l=1:NP
    disp(['mapping vector field..', num2str(100*l/NP), '%']);
    for m=1:NP
        for n=1:NP

            inv_y = locally_weighted_translation_inverse(rho,c,v,[y1(l,m,n);y2(l,m,n);y3(l,m,n)]);
            J = locally_weighted_translation_derivative(rho,c,v,inv_y);
            ydot = J*A*(inv_y - qf);
            
            y1dot(l,m,n) = ydot(1);
            y2dot(l,m,n) = ydot(2);
            y3dot(l,m,n) = ydot(3);
        end
    end
end

subplot(1,2,2); hold on;
axis equal; axis([0 1 0 1]);
axis equal; axis([plot_min(1) plot_max(1) plot_min(2) plot_max(2) plot_min(3) plot_max(3)]);
view([135 35]); xlabel('x'); ylabel('y'); zlabel('z');

streamline(x1,x2,x3,y1dot,y2dot,y3dot,startpoints(1,:),startpoints(2,:),startpoints(3,:));

%%
task.Xd = Y;
task.rho = rho;
task.c = c;
task.v = v;
task.A = A;
task.qi = qi;
task.qf = qf;
task.T_d = T_d;

%%
disp('done');