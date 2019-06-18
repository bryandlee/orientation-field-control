close all
clear
clc
%% Settings
t = 0;
dt = 1e-2;
horizon = 10;
export_video = true;

%% Velocity Field
load('../data/vf1.mat');

% NP = 100;    % num points for vectorfield plot
% vf.x1 = linspace(-1,1,NP);
% vf.x2 = vf.x1;
% vf.x3 = [0 1];
% [vf.x1, vf.x2, vf.x3] = meshgrid(vf.x1, vf.x2, vf.x3);
% 
% eigen_direction = (vf.qf-vf.qi)/norm(vf.qf-vf.qi);
% normal_direction = [eigen_direction(2); -eigen_direction(1)];
% vf.A = [eigen_direction normal_direction]' * [-1 0; 0 -5] * [eigen_direction normal_direction];
% 
% x1dot = vf.A(1,1)*(vf.x1 - vf.qf(1)) + vf.A(1,2)*(vf.x2 - vf.qf(2));
% x2dot = vf.A(2,1)*(vf.x1 - vf.qf(1)) + vf.A(2,2)*(vf.x2 - vf.qf(2));
% x3dot = - (vf.x3 - 0);
% 
% y1 = vf.x1;
% y2 = vf.x2;
% y3 = vf.x3;
% vf.x1dot = x1dot;
% vf.x2dot = x2dot;
% vf.x3dot = x3dot;
% 
% for l=1:NP
%     disp(['mapping vector field..', num2str(100*l/NP), '%']);
%     for m=1:NP
%         for n=1:1
% 
%             inv_y = locally_weighted_translation_inverse(vf.rho,vf.c,vf.v,[y1(l,m,n);y2(l,m,n)]);
%             J = locally_weighted_translation_derivative(vf.rho,vf.c,vf.v,inv_y);
%             ydot = J*vf.A*(inv_y - vf.qf);
%             
%             vf.x1dot(l,m,n) = ydot(1);
%             vf.x2dot(l,m,n) = ydot(2);
%         end
%     end
% end

figure(1); hold on;
axis equal; axis([-1 1 -1 1]);
plot(vf.Xd(1,:), vf.Xd(2,:),'-','Color',[0.5,0.5,0.5],'LineWidth',2);
% quiver(x1, x2, y1dot, y2dot ,'r');
streamplot = streamslice(vf.x1,vf.x2,vf.x3,vf.x1dot,vf.x2dot,vf.x3dot,[],[],0);
set(streamplot,'Color',[0.5 0.5 0.5]);

%% Video export
    if export_video
        writerObj = VideoWriter('../data/example','MPEG-4'); % 
        writerObj.FrameRate = 15;
    
        % open the video writer
        open(writerObj);
    end

%% Init robot
robot = makePlanar3R();

dof = robot.dof;
qi = rand(dof,1)*3;

% plot
figure(1); hold on;
axis equal;
axis([-1 1 -1 1]);

% forward kinematics
T = zeros(4,4,dof+1);
for i = 1:robot.dof
    T(:,:,i) = solveForwardKinematics(qi(1:i), robot.A(:,1:i), robot.M(:,:,1:i));
end
T(:,:,dof+1) = T(:,:,dof)*robot.M_ee;
xdata = reshape(T(1,4,:),[1,dof+1]);
ydata = reshape(T(2,4,:),[1,dof+1]);

links = plot(xdata, ydata, 'b', 'LineWidth', 2);
joints = plot(xdata(1:end-1), ydata(1:end-1), '.k', 'MarkerSize', 10);

J_b = getJacobianBody(zeros(dof,1), robot.A_b);
V_b = J_b*zeros(dof,1);
V_sb = T(1:3,1:3,end)*V_b(4:6);
end_effector_V = plot([xdata(end) xdata(end)+V_sb(1)], [ydata(end) ydata(end)+V_sb(2)], 'r', 'LineWidth', 1);

desired_V_sb = V_sb;
end_effector_V_desired = plot([xdata(end) xdata(end)+desired_V_sb(1)], [ydata(end) ydata(end)+desired_V_sb(2)], 'g', 'LineWidth', 1);

getframe;

%% Simulation
q = qi;
qdot = zeros(dof,1);
tau = zeros(dof,1);
timer = text(1,1,num2str(t));

while t < horizon
   
    %%%%%%%%%%% Control %%%%%%%%%%%%
    
    X_sb = T(1:2,4,end);
    inv_X_sb = locally_weighted_translation_inverse(vf.rho,vf.c,vf.v,X_sb);
    J = locally_weighted_translation_derivative(vf.rho,vf.c,vf.v,inv_X_sb);

    desired_V_sb(1:2) = J*vf.A*(inv_X_sb - vf.qf);
    
    desired_V_b = [0;0;0 ; T(1:3,1:3,end)'*desired_V_sb];

    
    lambda1 = 1;
    lambda2 = 1;
    
    e1 = normc(desired_V_b);
    Q = [e1, null(e1(:).')];
    Lambda = diag([0, lambda1*ones(1,6-1)]);
   
    D= Q*Lambda*Q';
    
    gravity_compensation = solveInverseDynamics(robot.A,robot.M,q,zeros(dof,1),zeros(dof,1),robot.G,[0;0;0;0;9.8;0]);
    
    
    tau = gravity_compensation + J_b'*(-D*V_b + lambda2*desired_V_b);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    % Runge-Kutta 4th order integration
    [~, q_and_qdot] = ode45(@(t,q_and_qdot) robot_dynamics(t,q_and_qdot,tau,robot), t:dt/2:t+dt, [q; qdot]);
    
    q = q_and_qdot(3,1:end/2)';
    qdot = q_and_qdot(3,end/2+1:end)';

    % plot
    T = zeros(4,4,dof+1);
    for i = 1:robot.dof
        T(:,:,i) = solveForwardKinematics(q(1:i), robot.A(:,1:i), robot.M(:,:,1:i));
    end
    T(:,:,dof+1) = T(:,:,dof)*robot.M_ee;

    xdata = reshape(T(1,4,:),[1,dof+1]);
    ydata = reshape(T(2,4,:),[1,dof+1]);
    
    set(links, 'Xdata', xdata, 'Ydata', ydata);
    set(joints, 'Xdata', xdata(1:end-1), 'Ydata', ydata(1:end-1));
   
    J_b = getJacobianBody(q, robot.A_b);
    V_b = J_b*qdot;
    V_sb = T(1:3,1:3,end)*V_b(4:6);
    set(end_effector_V, 'Xdata', [xdata(end) xdata(end)+1e0*V_sb(1)], 'Ydata', [ydata(end) ydata(end)+1e0*V_sb(2)]);
    set(end_effector_V_desired, 'Xdata', [xdata(end) xdata(end)+1e0*desired_V_sb(1)], 'Ydata', [ydata(end) ydata(end)+1e0*desired_V_sb(2)]);

    set(timer, 'String', num2str(t));
    
    frame = getframe(gcf);
        
    if export_video
        writeVideo(writerObj, frame);
    end
    
    t = t + dt;
end

if export_video
    % close the writer object
    close(writerObj);
end

%% Functions

% ROBOT DYNAMICS - state equation
% q_and_qdot = [q; qdot] in R^2n, tau in R^n
function dXdt = robot_dynamics(t,q_and_qdot,tau,robot)
    q = q_and_qdot(1:end/2);
    qdot = q_and_qdot(end/2+1:end);
    
    c1 = 0.2; c2 = 0.2;
    friction = -c1*qdot - c2*sign(qdot);
    
    tau = tau + friction;
    
    dXdt = [qdot; solveForwardDynamics(robot.A,robot.M,q,qdot,tau,robot.G,[0;0;0;0;9.8;0])];
end

% 