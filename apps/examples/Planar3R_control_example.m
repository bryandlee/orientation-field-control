close all
clear
clc
%% Settings
t = 0;
dt = 1e-2;

%% Init robot
robot = makePlanar3R();

dof = robot.dof;
qi = rand(dof,1)*3;

% plot
figure(); hold on;
axis equal;
axis([-1 1 -1 1]);

T = zeros(4,4,dof+1);
for i = 1:robot.dof
    T(:,:,i) = solveForwardKinematics(qi(1:i), robot.A(:,1:i), robot.M(:,:,1:i));
end
T(:,:,dof+1) = T(:,:,dof)*robot.M_ee;
xdata = reshape(T(1,4,:),[1,dof+1]);
ydata = reshape(T(2,4,:),[1,dof+1]);

links = plot(xdata, ydata, 'b', 'LineWidth', 2);
joints = plot(xdata(1:end-1), ydata(1:end-1), '.k', 'MarkerSize', 10);

V_b = getJacobianBody(zeros(dof,1), robot.A_b)*zeros(dof,1);
V_sb = T(1:3,1:3,end)*V_b(4:6);
end_effector_V = plot([xdata(end) xdata(end)+V_sb(1)], [ydata(end) ydata(end)+V_sb(2)], 'r', 'LineWidth', 2);

getframe;

%% Simulation
q = qi;
qdot = zeros(dof,1);
tau = zeros(dof,1);
timer = text(1,1,num2str(t));

while true
    
    %%%%%%%%%%% Control %%%%%%%%%%%%
    
%     gravity_compensation = solveInverseDynamics(robot.A,robot.M,q,zeros(dof,1),zeros(dof,1),robot.G,[0;0;0;0;9.8;0]);
%     
%     
%     tau = gravity_compensation;    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Runge-Kutta 4th order integration
    [~, X] = ode45(@(t,X) robot_dynamics(t,X,tau,robot), t:dt/2:t+dt, [q; qdot]);

    q = X(3,1:end/2)';
    qdot = X(3,end/2+1:end)';

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
   
    V_b = getJacobianBody(q, robot.A_b)*qdot;
    V_sb = T(1:3,1:3,end)*V_b(4:6);
    set(end_effector_V, 'Xdata', [xdata(end) xdata(end)+1e-5*V_sb(1)], 'Ydata', [ydata(end) ydata(end)+1e-5*V_sb(2)]);

    set(timer, 'String', num2str(t));
    getframe;
    
    t = t + dt;
end

%% Functions

% ROBOT DYNAMICS - state equation
% X = [q; qdot] in R^2n, tau in R^n
function dXdt = robot_dynamics(t,X,tau,robot)
    q = X(1:end/2);
    qdot = X(end/2+1:end);
    
    c1 = 0.2; c2 = 0.2;
    friction = -c1*qdot - c2*sign(qdot);
    
    tau = tau + friction;
    
    dXdt = [qdot; solveForwardDynamics(robot.A,robot.M,q,qdot,tau,robot.G,[0;0;0;0;9.8;0])];
end

% 