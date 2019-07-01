close all
clear
clc

%% Settings
t = 0;
dt = 1e-2;
horizon = 10;
export_video = true;
gravity = [0;0;0;0;0;9.81];

% controller   
E = 50;
M_f = 100;
gamma = 0.5;
w_n = 5;
delta = 1;
eta = 1;

% perturbation
F1 = rand(6,1)*10;
F2 = rand(6,1);

%% Velocity Field
load('../data/task3.mat');

%% Init robot
robot = makeFrankaPanda();

dof = robot.dof;
qi = getFeasiblePose(robot);
tau_min = robot.tau_min;
tau_max = robot.tau_max;
dtau_max = robot.dtau_max;

%% 

% %% Desried Path Plot
% visualizer2 = visualizePandaRealtime(robot,qi);
% 
% visualizer2.fig;
% plot3(task.qf(1),task.qf(2),task.qf(3),'*b');
% for i = 1:100
%     plot_SE3(task.T_d(:,:,i));
% end
% 
% %% Integral Curve
% visualizer3 = visualizePandaRealtime(robot,qi);
% 
% visualizer3.fig;
% plot3(task.qf(1),task.qf(2),task.qf(3),'*b');
% for i = 1:100
%     plot_SE3(task.T_d(:,:,i));
% end
% L=1000;
% 
% T_init = solveForwardKinematics(qi, robot.A, robot.M)* robot.M_ee;
% T_e=zeros(4,4,L);
% T_e(:,:,1) = T_init;
% 
% q = qi;
% E_max = 0;
% 
% for i=1:L-1
%     T_current = T_e(:,:,i);
%     x_current = T_current(1:3,4);
%     R_current = T_current(1:3,1:3);
% 
%     T_current_ = inverse_SE3(task.T_zero)*T_current;
%     x_current_ = T_current_(1:3,4);
%     R_current_ = T_current_(1:3,1:3);
%     
%     R_zero= task.T_zero(1:3,1:3);
%     
%     v_bar = getVelocityField(x_current, task);
%     v_bar_ = R_zero'*v_bar;
% 
%     % [W_d,V_d,W_d_dot,V_d_dot] = PVFC_Velocity(x_current_,x_dot_current_,R_current_,R_dot_current_,Theta,index,v_bar_,v_bar_dot_,w_n);
%     R_g = OrientationFieldGen3D(x_current_,task.Theta,task.index);
%     % plot_SE3(T_zero*[R_g x_current_; 0 0 0 1])
% 
%     KK = K(R_g,R_current_,delta);
%     V_d = v_d(KK,v_bar_);
%     k = k_(R_g,R_current_,w_n,eta);
%     dR_g_v_d = dR_g(x_current_,V_d,task.Theta,task.index);
%     W_d = w_d(R_g,R_current_,k,dR_g_v_d);
% 
%     V_b = [R_current_'*skew(W_d) ; R_current_'*V_d];
%     T_e(:,:,i+1)=T_e(:,:,i)*exp_se3(V_b*dt);
% 
%     plot_SE3(T_e(:,:,i));
%     getframe;
%     i
% end
% 

%% 
visualizer = visualizePandaRealtime(robot,qi);

visualizer.fig;
plot3(task.qf(1),task.qf(2),task.qf(3),'*b');
for i = 1:size(task.T_d,3)
    plot_SE3(task.T_d(:,:,i));
end
attractor = text(task.qf(1)+0.05,task.qf(2)+0.05,task.qf(3)+0.05,' ');

V_sb_desired = zeros(6,1);
Vdot_sb_desired = zeros(6,1);
V_sb = zeros(6,1);
F_perturb = zeros(6,1);

T = solveForwardKinematics(qi, robot.A, robot.M) * robot.M_ee;
V_plot = plot3([T(1,4) T(1,4)+V_sb(4)], [T(2,4) T(2,4)+V_sb(5)], [T(3,4) T(3,4)+V_sb(6)], 'r', 'LineWidth', 1);
V_des_plot = plot3([T(1,4) T(1,4)+V_sb_desired(4)], [T(2,4) T(2,4)+V_sb_desired(5)], [T(3,4) T(3,4)+V_sb_desired(6)], 'g', 'LineWidth', 1);
F_perturb_plot = plot3([T(1,4) T(1,4)+F_perturb(4)], [T(2,4) T(2,4)+F_perturb(5)], [T(3,4) T(3,4)+F_perturb(6)], 'b', 'LineWidth', 2);


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
tau_prev = tau;
timer = text(1,1,1,num2str(t));

while t < horizon
       
    % Runge-Kutta 4th order integration
%     [~, q_and_qdot] = ode45(@(t,q_and_qdot) robot_dynamics(t,q_and_qdot,tau,robot), t:dt/2:t+dt, [q; qdot]);
    q_and_qdot = RK4(@(t,q_and_qdot) robot_dynamics(t,q_and_qdot,tau,robot), t, t+dt, dt/100, [q; qdot]);
    
    q = q_and_qdot(1:end/2);
    qdot = q_and_qdot(end/2+1:end);

    t = t + dt;

    %%%%%%%%%%% Control %%%%%%%%%%%%
    T = solveForwardKinematics(q, robot.A, robot.M) * robot.M_ee;
    
    J_b = getJacobianBody(q, robot.A_b);
    V_b = J_b * qdot;

    Ad_analytic = zeros(6,6);
    Ad_analytic(1:3,1:3) = T(1:3,1:3);
    Ad_analytic(4:6,4:6) = T(1:3,1:3);
    J_a = Ad_analytic*J_b;
    V_sb = J_a*qdot;
    
    % desired v and vdot
    X_sb = T(1:3,4);
    [v_sb_desired, vdot_sb_desired] = getVelocityFieldDerivative(X_sb, V_sb(4:6), task);

    % orientation field
    T_0b = inverse_SE3(task.T_zero)*T;
    x_0b = T_0b(1:3,4);
    R_0b = T_0b(1:3,1:3);
    
    xdot_0b = R_0b * V_b(4:6);
    Rdot_0b = R_0b * skew(V_b(1:3));
    
    R_0s = task.T_zero(1:3,1:3)';
    v_desired_0 = R_0s * v_sb_desired;
    vdot_desired_0 = R_0s * vdot_sb_desired;
    
    [W_0d,V_0d,W_0d_dot,V_0d_dot] = PVFC_Velocity(x_0b,xdot_0b,R_0b,Rdot_0b,task.Theta,task.index,v_desired_0,vdot_desired_0,w_n,delta,eta);
    
    R_s0 = R_0s';
    V_sb_desired(1:3) = R_s0 * skew(W_0d);
    V_sb_desired(4:6) = R_s0 * V_0d;
    Vdot_sb_desired(1:3) = R_s0 * skew(W_0d_dot);
    Vdot_sb_desired(4:6) = R_s0 * V_0d_dot;
    
    % Analytic Jacobian derivative
    dAd_analytic = zeros(6,6);
    dAd_analytic(1:3,1:3) = skew(V_b(1:3));
    dAd_analytic(4:6,4:6) = skew(V_b(1:3));
    dAd_analytic = Ad_analytic*dAd_analytic;
    dJ_b = getJacobianBodyDerivative(q, robot.A_b, qdot);
    dJ_a = dAd_analytic*J_b + Ad_analytic*dJ_b;
    
    % Jacobian inverse
    JJT_inv = pinv(J_a*J_a');
    J_inv = J_a'*JJT_inv;
    
    dJ_a_inv = - J_a' * JJT_inv*(dJ_a*J_a' + J_a*dJ_a')*JJT_inv + dJ_a'*JJT_inv;
    
    % v and vdot
    v_q = J_inv*V_sb_desired;
    vdot_q = dJ_a_inv* V_sb_desired + J_inv * Vdot_sb_desired;


    tau = getPVFC(robot.A,robot.M,q,qdot,v_q,vdot_q,robot.G,gravity,E,M_f,gamma);
    tau = tau(1:dof);
    
    % joint input torque constraint
    tau_max_cutoff = max(tau_max, tau_prev + dtau_max*dt);
    tau_max_indicies = tau > tau_max_cutoff;
    tau(tau_max_indicies) = tau_max_cutoff(tau_max_indicies);
    
    tau_min_cutoff = min(tau_min, tau_prev - dtau_max*dt);
    tau_min_indicies = tau < tau_min_cutoff;
    tau(tau_min_indicies) = tau_min_cutoff(tau_min_indicies);

    if (~all(~tau_max_indicies)) || (~all(~tau_min_indicies))
        disp(' ************* torque limit *************');
    end
    
    % perturbation
    if (1.0 < t) && (t < 1.1)
        F_perturb = F1;
    elseif (3.0 < t) && (t < 3.1)
        F_perturb = F2;
    else
        F_perturb = zeros(6,1);
    end
    
    tau = tau + J_b'*Ad_analytic'*F_perturb;

    
    % check status
    M_q = getMassMatrix(robot.A, robot.M, q, robot.G);
    kinetic_energy = qdot'*M_q*qdot/2
    kinetic_energy_desired = v_q'*M_q*v_q/2
    norm(v_q)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % plot
    set(timer, 'String', num2str(t));
    set(V_plot, 'Xdata', [T(1,4) T(1,4)+V_sb(4)], 'Ydata', [T(2,4) T(2,4)+V_sb(5)], 'Zdata', [T(3,4) T(3,4)+V_sb(6)]);
    set(V_des_plot, 'Xdata', [T(1,4) T(1,4)+V_sb_desired(4)], 'Ydata', [T(2,4) T(2,4)+V_sb_desired(5)], 'Zdata', [T(3,4) T(3,4)+V_sb_desired(6)]);
    set(F_perturb_plot, 'Xdata', [T(1,4) T(1,4)+F_perturb(4)*0.1], 'Ydata', [T(2,4) T(2,4)+F_perturb(5)*0.1], 'Zdata', [T(3,4) T(3,4)+F_perturb(6)*0.1]);
%     plot_SE3(T);

    visualizePandaRealtime(robot,q,visualizer);

    frame = getframe(gcf);
        
    if export_video
        writeVideo(writerObj, frame);
    end
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

