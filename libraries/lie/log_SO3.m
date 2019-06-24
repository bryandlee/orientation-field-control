function w = log_SO3( R )
    %LOG_SO3 이 함수의 요약 설명 위치
    %   자세한 설명 위치
    s = size(R);
    dim = s(1); 
    R=reshape(R,[dim,dim]);
    
    EPS = 1e-15;
    if dim == 3
        cos = (R(1,1) + R(2,2) + R(3,3) - 1)/2;
        w = zeros(3,1);
       
        if cos == -1
            if R(1,1) > 1 - EPS
                w(1) = pi;
                return
            elseif R(2,2) > 1 - EPS
                w(2) = pi;
                return
            elseif R(3,3) > 1 - EPS
                w(3) = pi;
                return
            end
            
            w(1) = sqrt((R(2,1) * R(2,1) + R(3,1) * R(3,1)) / (1.0 - R(1,1)));
            w(2) = sqrt((R(1,2) * R(1,2) + R(3,2) * R(3,2)) / (1.0 - R(2,2)));
            w(3) = sqrt((R(1,3) * R(1,3) + R(2,3) * R(2,3)) / (1.0 - R(3,3)));
            w = w * pi/sqrt(2);
            return
        end
        
        if cos > 1 
            cos = 1;
        elseif cos < -1
            cos = -1;
        end
        
        theta = acos(cos);
        if theta < EPS && theta > -EPS
            t_st = 3/(6 - theta*theta);
        else
            t_st = theta / (2*sin(theta));
        end
        
        w(1) = R(3,2)-R(2,3);
        w(2) = R(1,3)-R(3,1);
        w(3) = R(2,1)-R(1,2);
        w = w*t_st;
        return

    elseif dim == 2
        a = acos(R(1,1));
        b = -a;
        if sin(a)-(R(2,1)) < 1.0e-10
            theta = a;
        elseif sin(b)-(R(2,1)) < 1.0e-10
            theta = b;
        end
        w = [0,-theta;theta,0];
    end
end

