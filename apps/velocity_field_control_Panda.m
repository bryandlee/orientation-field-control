close all
clear
clc
%% Settings
t = 0;
dt = 1e-2;
horizon = 10;
export_video = false;

%% Velocity Field
load('../data/task3.mat');

%% Init robot
robot = makeFrankaPanda();

dof = robot.dof;
qi = getFeasiblePose(robot);
visualizer = visualizePandaRealtime(robot,qi);

visualizer.fig;
plot3(task.Xd(1,:), task.Xd(2,:), task.Xd(3,:),'-','Color',[0.5,0.5,0.5],'LineWidth',2);
plot3(task.qf(1),task.qf(2),task.qf(3),'*b');
attractor = text(task.qf(1)+0.05,task.qf(2)+0.05,task.qf(3)+0.05,'attractor');

desired_V_sb = zeros(6,1);

%% Video export
if export_video
    writerObj = VideoWriter('../data/example_panda','MPEG-4'); % 
    writerObj.FrameRate = 15;

    % open the video writer
    open(writerObj);
end

%% Simulation
q = qi;
qdot = zeros(dof,1);
tau = zeros(dof,1);
timer = text(1,1,1,num2str(t));

while t < horizon
   
    %%%%%%%%%%% Control %%%%%%%%%%%%
    T = solveForwardKinematics(q, robot.A, robot.M) * robot.M_ee;
    
    J_b = getJacobianBody(q, robot.A_b);
    V_b = J_b*qdot;
    Ad_analytic = zeros(6,6);
    Ad_analytic(1:3,1:3) = T(1:3,1:3);
    Ad_analytic(4:6,4:6) = T(1:3,1:3);
    V_sb = Ad_analytic*V_b;
    
    X_sb = T(1:3,4);
    inv_X_sb = locally_weighted_translation_inverse(task.rho,task.c,task.v,X_sb);
    J = locally_weighted_translation_derivative(task.rho,task.c,task.v,inv_X_sb);

    desired_V_sb(4:6) = J*task.A*(inv_X_sb - task.qf);
    desired_V_b = Ad_analytic' * desired_V_sb;
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    % Runge-Kutta 4th order integration
%     [~, q_and_qdot] = ode45(@(t,q_and_qdot) robot_dynamics(t,q_and_qdot,tau,robot), t:dt/2:t+dt, [q; qdot]);
    q_and_qdot = RK4(@(t,q_and_qdot) robot_dynamics(t,q_and_qdot,tau,robot), t, t+dt, dt/100, [q; qdot]);
    
    q = q_and_qdot(1:end/2);
    qdot = q_and_qdot(end/2+1:end);

    % plot
    visualizePandaRealtime(robot,q,visualizer);
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
    
    dXdt = [qdot; solveForwardDynamics(robot.A,robot.M,q,qdot,tau,robot.G,[0;0;0;0;0;9.81])];
end

% Runge Kutta Method 4th Order 
function y= RK4(f,t0,t1,dt,y0)
    t = t0:dt:t1;
    y = y0;

    for i=1:(length(t)-1)
        k1 = f(t(i),          y); 
        k2 = f(t(i) + 0.5*dt, y + 0.5*dt*k1); 
        k3 = f(t(i) + 0.5*dt, y + 0.5*dt*k2); 
        k4 = f(t(i) + dt,     y + k3*dt); 

        y = y + (1/6)*(k1+2*k2+2*k3+k4)*dt;
    end
end

