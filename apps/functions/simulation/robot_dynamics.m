% ROBOT DYNAMICS - state equation
% q_and_qdot = [q; qdot] in R^2n, tau in R^n
function dXdt = robot_dynamics(t,q_and_qdot,tau,robot,gravity)
    q = q_and_qdot(1:end/2);
    qdot = q_and_qdot(end/2+1:end);
    
    c1 = 0.2; c2 = 0.2;
    friction = -c1*qdot - c2*sign(qdot);
    
    tau = tau + friction;
    
    dXdt = [qdot; solveForwardDynamics(robot.A,robot.M,q,qdot,tau,robot.G,gravity)];
end

