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