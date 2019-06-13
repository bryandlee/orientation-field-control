robot = makePlanar3R();

dof = robot.dof;

% random spline trajectory
complexity  = 3;   % num of spline coefs per joints
basis_order = 4;   % cubic spline
horizon  = 20;      % trajectory horizon

num_samples = 500;
spline_params = rand(complexity,dof);
qi = zeros(dof,1);
qf = rand(dof,1)*10;
t = linspace(0,horizon,num_samples);
[q, qdot] = makeSplineP2P(qi,qf,spline_params, basis_order, horizon, t);


% plot
figure(); hold on;
axis equal;
axis([-1 1 -1 1]);

T = zeros(4,4,dof+1);
for i = 1:robot.dof
    T(:,:,i) = solveForwardKinematics(zeros(i,1), robot.A(:,1:i), robot.M(:,:,1:i));
end
T(:,:,dof+1) = T(:,:,dof)*robot.M_ee;
x = reshape(T(1,4,:),[1,dof+1]);
y = reshape(T(2,4,:),[1,dof+1]);

links = plot(x, y, 'b', 'LineWidth', 2);
joints = plot(x(1:end-1), y(1:end-1), '.k', 'MarkerSize', 10);
end_effector_T = plot_SE3(T(:,:,end));

V_b = getJacobianBody(zeros(dof,1), robot.A_b)*zeros(dof,1);
V_bs = T(1:3,1:3,end)*V_b(4:6);
end_effector_V = plot([x(end) x(end)+V_bs(1)], [y(end) y(end)+V_bs(2)], 'r', 'LineWidth', 2);

getframe;
    
for t = 1:num_samples
    T = zeros(4,4,dof+1);
    for i = 1:robot.dof
        T(:,:,i) = solveForwardKinematics(q(1:i,t), robot.A(:,1:i), robot.M(:,:,1:i));
    end
    T(:,:,dof+1) = T(:,:,dof)*robot.M_ee;

    x = reshape(T(1,4,:),[1,dof+1]);
    y = reshape(T(2,4,:),[1,dof+1]);
    
    set(links, 'Xdata', x, 'Ydata', y);
    set(joints, 'Xdata', x(1:end-1), 'Ydata', y(1:end-1));
    plot_SE3(T(:,:,end), end_effector_T);
    
    V_b = getJacobianBody(q(:,t), robot.A_b)*qdot(:,t);
    V_bs = T(1:3,1:3,end)*V_b(4:6);
    set(end_effector_V, 'Xdata', [x(end) x(end)+V_bs(1)], 'Ydata', [y(end) y(end)+V_bs(2)]);

    getframe;

end