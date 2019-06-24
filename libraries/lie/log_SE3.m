function S = log_SE3( T )
    R = T(1:3,1:3);
    p = T(1:3,4);

    S = zeros(6,1);
    
    EPS = 1e-12;
    cost = (R(1,1) + R(2,2) + R(3,3) - 1)/2;
    w = zeros(3,1);

    if cost == -1
        if R(1,1) > 1 - EPS
            w(1) = pi;
        elseif R(2,2) > 1 - EPS
            w(2) = pi;
        elseif R(3,3) > 1 - EPS
            w(3) = pi;
        else
            w(1) = sqrt((R(2,1) * R(2,1) + R(3,1) * R(3,1)) / (1.0 - R(1,1)));
            w(2) = sqrt((R(1,2) * R(1,2) + R(3,2) * R(3,2)) / (1.0 - R(2,2)));
            w(3) = sqrt((R(1,3) * R(1,3) + R(2,3) * R(2,3)) / (1.0 - R(3,3)));
            w = w * pi/sqrt(2);         
        end
        
        w2 = zeros(3,1);
        w2(1) = w(1)*w(1);
        w2(2) = w(2)*w(2);
        w2(3) = w(3)*w(3);

        w3 = zeros(3,1);
        w3(1) = w2(1)*w(1);
        w3(2) = w2(2)*w(2);
        w3(3) = w2(3)*w(3);
        
        v = zeros(3,1);
        v(1) = w(2)*w(3);
        v(2) = w(3)*w(1);
        v(3) = w(1)*w(2);

        id = 0.25 * pi * pi / (w2(1)*w2(1) + w2(2)*w2(2) + w2(3)*w2(3) +  2*(w2(1)*w2(2) + w2(2)*w2(3) + w2(3)*w2(1)));
        p = p * id;
        
        S(1:3) = w;
        S(4) = 2 * (2 * w2(1) * p(1) + (w3(3) + v(2) * w(1) + v(1) * w(2) + 2 * v(3)) * p(2) + (2 * v(2) - w3(2) - w(1) * v(3) - w(3) * v(1)) * p(3));
        S(5) = 2 * ((2 * v(3) - w3(3) - v(1) * w(2) - v(2) * w(1)) * p(1) + 2 * w2(2) * p(2) + (w3(1) + w(3) * v(2) + v(3) * w(2) + 2 * v(1)) * p(3));
        S(6) = 2 * ((w(1) * v(3) + w(3) * v(1) + 2 * v(2) + w3(2)) * p(1) + (2 * v(1) - w3(1) - v(3) * w(2) - w(3) * v(2)) * p(2) + 2 * w2(3) * p(3));
        
        return
    else
        if cost > 1 
            cost = 1;
        elseif cost < -1
            cost = -1;
        end
        
        theta = acos(cost);
        if theta < EPS && theta > -EPS
            t_st = 3/(6 - theta*theta);
        else
            t_st = theta / (2*sin(theta));
        end
        
        w(1) = R(3,2)-R(2,3);
        w(2) = R(1,3)-R(3,1);
        w(3) = R(2,1)-R(1,2);

        w2 = zeros(3,1);
        w2(1) = w(1)*w(1);
        w2(2) = w(2)*w(2);
        w2(3) = w(3)*w(3);

        w3 = zeros(3,1);
        w3(1) = w(2)*w(3);
        w3(2) = w(3)*w(1);
        w3(3) = w(1)*w(2);

        w = w*t_st;
        
        st = sin(theta);

        if theta < EPS && theta > -EPS
            stct = theta / 48;
        else
            stct = (2*st - theta*(1+cos(theta)))/(8*st*st*st);
        end        

        S(1:3) = w;
        S(4) = (1 - stct * (w2(2) + w2(3))) * p(1) + (0.5 * w(3) + stct * w3(3)) * p(2) + (stct * w3(2) - 0.5 * w(2)) * p(3);
        S(5) = (stct * w3(3) - 0.5 * w(3)) * p(1) + (1 - stct * (w2(1) + w2(3))) * p(2) + (0.5 * w(1) + stct * w3(1)) * p(3);
        S(6) = (0.5 * w(2) + stct * w3(2)) * p(1) + (stct * w3(1) - 0.5 * w(1)) * p(2) + (1 - stct * (w2(2) + w2(1))) * p(3);        
        return
    end
    
end
