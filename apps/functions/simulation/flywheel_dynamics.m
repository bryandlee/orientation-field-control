% Virtual Flywheel Dynamics
function qddot_flywheel = flywheel_dynamics(t,tau,M_f)
    qddot_flywheel = tau/M_f;
end