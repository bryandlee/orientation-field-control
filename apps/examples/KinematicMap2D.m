close all
clear 
clc

%% Random Path on SE(2) Generation
%% Settings
dim = 2;    % data dimension
N = 50;    % num of data

K = 150;    % num of transformations
mu = 0.9;   % safety margin
beta = 0.5; % deformation rate

%% Data generation
% data settings
complexity  = 1;   % num of spline coefs per joints
complexity_ori = 3;
basis_order = 4;   % cubic spline
horizon  = 1;      % trajectory horizon

qi = 5*rand(dim,1);  % init point
qf = 5*rand(dim,1);  % final point
qi_ori = 5*rand(dim*(dim-1)/2,1);  % init point
qf_ori = 5*rand(dim*(dim-1)/2,1);  % final point

spline_params = rand(complexity,dim);
spline_params_ori = rand(complexity_ori,dim*(dim-1)/2);

t = linspace(0,horizon,1e3*N);
[sp, dsp] = makeSplineP2P(qi,qf,spline_params, basis_order, horizon,t);
spline_velocity = vecnorm(dsp);
total_length = sum(spline_velocity);

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
% s = linspace(0,horizon,N);
Y = makeSplineP2P(qi,qf,spline_params, basis_order, horizon, s);
so3 = makeSplineP2P(qi_ori,qf_ori,spline_params_ori, basis_order, horizon, s);


% trajectory plot
figure('Name','Sample Trajectory','NumberTitle','off','units','pixels','pos',[100 100 800 800]);
subplot(2,2,1);
hold on; axis equal; 
title('Desired path in SE(2)')
% view([-135 25]);

if dim == 3
    plot3(X(1,:), X(2,:), X(3,:), '-k');
    plot3(Y(1,:), Y(2,:), Y(3,:), '--k');
    for i = 1:N
    T = [exp_so3(so3(:,i)) Y(:,i);
         0 0 0 1];
    plot_SE2(T);
    end
    
elseif dim == 2
    T_traj = zeros(3,3,N);
    T_zero = [exp_so2(so3(:,1)) Y(:,1);
         0 0 1]; 
    for i = 1:N
    T = [exp_so2(so3(:,i)) Y(:,i);
         0 0 1];
    T_traj(:,:,i)=inv(T_zero)*T;
    plot_SE2(T_traj(:,:,i));
    end    
    plot(reshape(T_traj(1,3,:),[N,1]), reshape(T_traj(2,3,:),[N,1]), '-k');
end

x2 = max(reshape(T_traj(1,3,:),[N,1]))+0.2;
x1 = min(reshape(T_traj(1,3,:),[N,1]))-0.2;
y2 = max(reshape(T_traj(2,3,:),[N,1]))+0.2;
y1 = min(reshape(T_traj(2,3,:),[N,1]))-0.2;
axis([x1 x2 y1 y2])

%% Connection Learning (W=0.00001)
n=5;
if dim ==3
    M = (n^3+6*n^2+11*n+6)/6;
elseif dim==2
    M = (n+2)*(n+1)/2;
end
W = 1e-5*eye(M);
%     k=2;
%     for i=1:n
%         for j=1:(i+1)*(i+2)/2
%             W(k,k) = i;
%             k=k+1;
%         end
%     end
[A,B,Theta] = Cal_ABTheta(T_traj,W,n);
R_g = zeros(dim,dim,N);
for i=1:N
    R_g(:,:,i) = OrientationFieldGen(reshape(T_traj(1:dim,dim+1,i),[dim,1]),Theta,n);
end
error = 0;
for i=1:N
    R_true = reshape(T_traj(1:dim,1:dim,i),[dim,dim]);
    error = error + norm(log_SO3(reshape(R_g(:,:,i),[dim,dim])' * R_true));
end
% err_avg(n+1) = error/N;

disp (error/N);

subplot(2,2,2);
hold on; axis equal; 
title('Orientation field (W=1e-5*Id)')
axis([x1 x2 y1 y2])
if dim ==3
    for i = 1:N
    T = [reshape(R_g(:,:,i),[dim,dim]) reshape(T_traj(1:dim,dim+1,i),[dim,1]);0 0 0 1];
    plot_SE2(T);
    end
elseif dim==2
    for i = 1:N
    T = [reshape(R_g(:,:,i),[dim,dim]) reshape(T_traj(1:dim,dim+1,i),[dim,1]);0 0 1];
    plot_SE2(T);
    end
    plot(reshape(T_traj(1,3,:),[N,1]), reshape(T_traj(2,3,:),[N,1]), '-k');
end

M = 5;
for i=1:N
    X_i  = reshape(T_traj(1:dim,dim+1,i),[dim,1]);
    if i < N
        v = reshape(T_traj(1:dim,dim+1,i+1),[dim,1])-reshape(T_traj(1:dim,dim+1,i),[dim,1]);
        V = (eye(dim)-v*v')*rand(dim,1);
        e = V/norm(V);
    elseif i == N
        v = reshape(T_traj(1:dim,dim+1,i),[dim,1])-reshape(T_traj(1:dim,dim+1,i-1),[dim,1]);
        V = (eye(dim)-v*v')*rand(dim,1);
        e = V/norm(V);
    end
    for j=1:M
    X = X_i + 0.3* (-1+2*j/(M+1))*e;
    R_G = OrientationFieldGen(X,Theta,n);
    if dim ==3
        T = [reshape(R_G,[dim,dim]) reshape(X,[dim,1]);0 0 0 1];
    elseif dim==2
        T = [reshape(R_G,[dim,dim]) reshape(X,[dim,1]);0 0 1];
    end
    plot_SE2(T);
    end
end
%% Connection Learning (W=0.001)
n=5;
if dim ==3
    M = (n^3+6*n^2+11*n+6)/6;
elseif dim==2
    M = (n+2)*(n+1)/2;
end
W = 0.001*eye(M);
%     k=2;
%     for i=1:n
%         for j=1:(i+1)*(i+2)/2
%             W(k,k) = i;
%             k=k+1;
%         end
%     end
[A,B,Theta] = Cal_ABTheta(T_traj,W,n);
R_g = zeros(dim,dim,N);
for i=1:N
    R_g(:,:,i) = OrientationFieldGen(reshape(T_traj(1:dim,dim+1,i),[dim,1]),Theta,n);
end
error = 0;
for i=1:N
    R_true = reshape(T_traj(1:dim,1:dim,i),[dim,dim]);
    error = error + norm(log_SO3(reshape(R_g(:,:,i),[dim,dim])' * R_true));
end
% err_avg(n+1) = error/N;

disp (error/N);

subplot(2,2,3);
hold on; axis equal; 
title('Orientation field (W=1e-3*Id)')
axis([x1 x2 y1 y2])
if dim ==3
    for i = 1:N
    T = [reshape(R_g(:,:,i),[dim,dim]) reshape(T_traj(1:dim,dim+1,i),[dim,1]);0 0 0 1];
    plot_SE2(T);
    end
elseif dim==2
    for i = 1:N
    T = [reshape(R_g(:,:,i),[dim,dim]) reshape(T_traj(1:dim,dim+1,i),[dim,1]);0 0 1];
    plot_SE2(T);
    end
    plot(reshape(T_traj(1,3,:),[N,1]), reshape(T_traj(2,3,:),[N,1]), '-k');
end

M = 5;
for i=1:N
    X_i  = reshape(T_traj(1:dim,dim+1,i),[dim,1]);
    if i < N
        v = reshape(T_traj(1:dim,dim+1,i+1),[dim,1])-reshape(T_traj(1:dim,dim+1,i),[dim,1]);
        V = (eye(dim)-v*v')*rand(dim,1);
        e = V/norm(V);
    elseif i == N
        v = reshape(T_traj(1:dim,dim+1,i),[dim,1])-reshape(T_traj(1:dim,dim+1,i-1),[dim,1]);
        V = (eye(dim)-v*v')*rand(dim,1);
        e = V/norm(V);
    end
    for j=1:M
    X = X_i + 0.3* (-1+2*j/(M+1))*e;
    R_G = OrientationFieldGen(X,Theta,n);
    if dim ==3
        T = [reshape(R_G,[dim,dim]) reshape(X,[dim,1]);0 0 0 1];
    elseif dim==2
        T = [reshape(R_G,[dim,dim]) reshape(X,[dim,1]);0 0 1];
    end
    plot_SE2(T);
    end
end
%% Connection Learning (W=0.1)
n=5;
if dim ==3
    M = (n^3+6*n^2+11*n+6)/6;
elseif dim==2
    M = (n+2)*(n+1)/2;
end
W = 0.1*eye(M);
%     k=2;
%     for i=1:n
%         for j=1:(i+1)*(i+2)/2
%             W(k,k) = i;
%             k=k+1;
%         end
%     end
[A,B,Theta] = Cal_ABTheta(T_traj,W,n);
R_g = zeros(dim,dim,N);
for i=1:N
    R_g(:,:,i) = OrientationFieldGen(reshape(T_traj(1:dim,dim+1,i),[dim,1]),Theta,n);
end
error = 0;
for i=1:N
    R_true = reshape(T_traj(1:dim,1:dim,i),[dim,dim]);
    error = error + norm(log_SO3(reshape(R_g(:,:,i),[dim,dim])' * R_true));
end
% err_avg(n+1) = error/N;

disp (error/N);


subplot(2,2,4);
hold on; axis equal; 
title('Orientation field (W=1e-1*Id)')
axis([x1 x2 y1 y2])
if dim ==3
    for i = 1:N
    T = [reshape(R_g(:,:,i),[dim,dim]) reshape(T_traj(1:dim,dim+1,i),[dim,1]);0 0 0 1];
    plot_SE2(T);
    end
elseif dim==2
    for i = 1:N
    T = [reshape(R_g(:,:,i),[dim,dim]) reshape(T_traj(1:dim,dim+1,i),[dim,1]);0 0 1];
    plot_SE2(T);
    end
    plot(reshape(T_traj(1,3,:),[N,1]), reshape(T_traj(2,3,:),[N,1]), '-k');
end

M = 5;
for i=1:N
    X_i  = reshape(T_traj(1:dim,dim+1,i),[dim,1]);
    if i < N
        v = reshape(T_traj(1:dim,dim+1,i+1),[dim,1])-reshape(T_traj(1:dim,dim+1,i),[dim,1]);
        V = (eye(dim)-v*v')*rand(dim,1);
        e = V/norm(V);
    elseif i == N
        v = reshape(T_traj(1:dim,dim+1,i),[dim,1])-reshape(T_traj(1:dim,dim+1,i-1),[dim,1]);
        V = (eye(dim)-v*v')*rand(dim,1);
        e = V/norm(V);
    end
    for j=1:M
    X = X_i + 0.3* (-1+2*j/(M+1))*e;
    R_G = OrientationFieldGen(X,Theta,n);
    if dim ==3
        T = [reshape(R_G,[dim,dim]) reshape(X,[dim,1]);0 0 0 1];
    elseif dim==2
        T = [reshape(R_G,[dim,dim]) reshape(X,[dim,1]);0 0 1];
    end
    plot_SE2(T);
    end
end

%% Veclocity filed on R^2




