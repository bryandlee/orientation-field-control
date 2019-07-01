function [v, vdot] = getVelocityFieldDerivative(x, xdot, task)

    dt = 1e-4;

    % v
    v = getVelocityField(x, task);

    % vdot (numerical)
    v_dv = getVelocityField(x + xdot*dt, task);
    v_dv_ = getVelocityField(x - xdot*dt, task);
    
    vdot = (v_dv - v_dv_)/(2*dt);
end