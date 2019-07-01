clear all;
clc;
close all;

load('../data/3.mat');
load('../data/task3.mat');
%% Settings
t = 0;
dt = 1e-2;
horizon = 10;
export_video = false;
%% Init robot
robot = makeFrankaPanda();

dof = robot.dof;
q_init = getFeasiblePose(robot);
qi = q_init;
visualizer = visualizePandaRealtime(robot,qi);

visualizer.fig;
plot3(task.Xd(1,:), task.Xd(2,:), task.Xd(3,:),'-','Color',[0.5,0.5,0.5],'LineWidth',2);
plot3(task.qf(1),task.qf(2),task.qf(3),'*b');
% attractor = text(task.qf(1)+0.05,task.qf(2)+0.05,task.qf(3)+0.05,'attractor');

%% Video export
if export_video
    writerObj = VideoWriter('../data/example_panda','MPEG-4'); % 
    writerObj.FrameRate = 15;

    % open the video writer
    open(writerObj);
end

%% Connection Learning
dim = 3;
poly_order = 3;
index = Indenx3DGen(poly_order);
ss = size(index);
N = ss(1);
W = 1.0e-1*eye(N);


T_zero = T_d(:,:,1);
R_zero = T_zero(1:3,1:3);
for i = 1:100
T_traj(:,:,i)=inv(T_zero)*T_d(:,:,i);
plot_SE3(T_d(:,:,i))
end
[A,B,Theta] = Cal_3DABTheta(T_traj,W,index);

R_g = zeros(dim,dim,100);

for i=1:100
    R_g(:,:,i) = OrientationFieldGen3D(reshape(T_traj(1:dim,dim+1,i),[dim,1]),Theta,index);
end
error = 0;
for i=1:100
    R_true = reshape(T_traj(1:dim,1:dim,i),[dim,dim]);
    error = error + norm(log_SO3(reshape(R_g(:,:,i),[dim,dim])' * R_true));
end
disp (error/100);

%% Orientation Field Visualization
% trajectory plot
figure('Name','Sample Trajectory','NumberTitle','off','units','pixels','pos',[100 100 800 800]);
hold on; axis equal; view([-135 25]);
N=100;
for i = 1:N
    X_i  = reshape(T_traj(1:dim,dim+1,i),[dim,1]);
    if i < N
        v = reshape(T_traj(1:dim,dim+1,i+1),[dim,1])-reshape(T_traj(1:dim,dim+1,i),[dim,1]);
        V = (eye(dim)-v*v')*rand(dim,1);
        e1 = V/norm(V);
        W = (eye(dim)-e1*e1')*(eye(dim)-v*v')*rand(dim,1);
        e2 = W/norm(W);
    elseif i == N
        v = reshape(T_traj(1:dim,dim+1,i),[dim,1])-reshape(T_traj(1:dim,dim+1,i-1),[dim,1]);
        V = (eye(dim)-v*v')*rand(dim,1);
        e1 = V/norm(V);
        W = (eye(dim)-e1*e1')*(eye(dim)-v*v')*rand(dim,1);
        e2 = W/norm(W);
    end
    M=3;
    for j=1:M
        for k=1:M
            X = X_i + 0.3* ((-1+2*j/(M+1))*e1 + (-1+2*k/(M+1))*e2);
            R_G = OrientationFieldGen3D(X,Theta,index);
            T = [reshape(R_G,[dim,dim]) reshape(X,[dim,1]);0 0 0 1];
            plot_SE3(T_zero*T);
        end
    end
end    

%% 
visualizer = visualizePandaRealtime(robot,qi);

visualizer.fig;
plot3(task.Xd(1,:), task.Xd(2,:), task.Xd(3,:),'-','Color',[0.5,0.5,0.5],'LineWidth',2);
plot3(task.qf(1),task.qf(2),task.qf(3),'*b');
% attractor = text(task.qf(1)+0.05,task.qf(2)+0.05,task.qf(3)+0.05,'attractor');
for i = 1:100
    T_traj(:,:,i)=inv(T_zero)*T_d(:,:,i);
    plot_SE3(T_d(:,:,i));
end

%% Integral Curve
w_n = 10;
delta = 1;
eta = 1;
L=1000;

T_init = solveForwardKinematics(q_init, robot.A, robot.M)* robot.M_ee;
T_e=zeros(4,4,L);
T_e(:,:,1) = T_init;

q = q_init;
E_max = 0;

for i=1:L-1
    T_current = T_e(:,:,i);
    x_current = T_current(1:3,4);
    R_current = T_current(1:3,1:3);

    T_current_ = inverse_SE3(T_zero)*T_current;
    x_current_ = T_current_(1:3,4);
    R_current_ = T_current_(1:3,1:3);

    v_bar = getVelocityField(x_current, task);
    v_bar_ = R_zero'*v_bar;

    % [W_d,V_d,W_d_dot,V_d_dot] = PVFC_Velocity(x_current_,x_dot_current_,R_current_,R_dot_current_,Theta,index,v_bar_,v_bar_dot_,w_n);
    R_g = OrientationFieldGen3D(x_current_,Theta,index);
    % plot_SE3(T_zero*[R_g x_current_; 0 0 0 1])

    KK = K(R_g,R_current_,delta);
    V_d = v_d(KK,v_bar_);
    k = k_(R_g,R_current_,w_n,eta);
    dR_g_v_d = dR_g(x_current_,V_d,Theta,index);
    W_d = w_d(R_g,R_current_,k,dR_g_v_d);

    V_b = [R_current_'*skew(W_d) ; R_current_'*V_d];
    T_e(:,:,i+1)=T_e(:,:,i)*exp_se3(V_b*dt);

    
    q = solveInverseKinematics(robot.A_s, T_e(:,:,i)*inverse_SE3(robot.M_sb), q);
    qdot_desired = pinv(getJacobianBody(q, robot.A_b)) * V_b;
    E = qdot_desired'* getMassMatrix(robot.A, robot.M, q, robot.G)* qdot_desired/2
    if E > E_max
        E_max = E;
    end
    
    visualizer.fig;
    plot_SE3(T_e(:,:,i));
    getframe;
end

