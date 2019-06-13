close all
clear 
clc

%% Settings
dim = 3;    % data dimension
N = 20;    % num of data

K = 150;    % num of transformations
mu = 0.9;   % safety margin
beta = 0.5; % deformation rate

%% Data generation
% data settings
complexity  = 3;   % num of spline coefs per joints
basis_order = 4;   % cubic spline
horizon  = 1;      % trajectory horizon

qi = 5*rand(dim,1);  % init point
qf = 5*rand(dim,1);  % final point
spline_params = rand(complexity,dim);
spline_params_ori = rand(complexity,dim);


t = linspace(0,horizon,1e3*N);
[sp, dsp] = makeSplineP2P(qi,qf,spline_params, basis_order, horizon,t);
spline_velocity = vecnorm(dsp);
total_length = sum(spline_velocity);

% constant velocity reparameterize
s = zeros(1,N); s(N) = horizon; j = 1;
for i = 2:N-1
    local_length = 0; k = 0;
    while true
        local_length = local_length + spline_velocity(j+k);
        if local_length > total_length/(N-1)
            s(i) = t(j+k-1);
            j = j+k-1;
            break;
        else
            k = k+1;
        end
    end
end

Y = makeSplineP2P(qi,qf,spline_params, basis_order, horizon, s);
so3 = makeSplineP2P(qi,qf,spline_params_ori, basis_order, horizon, s);

X = zeros(dim,N);
for i = 1:dim
    X(i,:) = linspace(Y(i,1),Y(i,N),N);
end

% trajectory plot
figure('Name','Sample Trajectory','NumberTitle','off','units','pixels','pos',[100 100 800 800]);
hold on; axis equal; view([-135 25]);
plot3(X(1,:), X(2,:), X(3,:), '-k');
plot3(Y(1,:), Y(2,:), Y(3,:), '--k');
for i = 1:N
    T = [exp_so3(so3(:,i)) Y(:,i);
         0 0 0 1];
    plot_SE3(T);
end
