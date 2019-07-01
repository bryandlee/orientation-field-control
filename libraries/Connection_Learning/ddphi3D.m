function DDPhi3D = ddphi3D(x,axis1,axis2,index)
x1 = x(1);
x2 = x(2);
x3 = x(3);
[N, b] = size(index);
DDPhi3D = zeros(N,1);
if axis1 == 1
    if axis2 == 1
        for i=1:N
        l = index(i,:);
        DDPhi3D(i,:) = l(1) * (l(1)-1) * x1^(l(1)-2) * x2^l(2) * x3^l(3);
        end
    elseif axis2 == 2
        for i=1:N
        l = index(i,:);
        DDPhi3D(i,:) = l(1) * l(2) * x1^(l(1)-1) * x2^(l(2)-1) * x3^l(3);
        end
    elseif axis2 == 3
        for i=1:N
        l = index(i,:);
        DDPhi3D(i,:) = l(1) * l(3) * x1^(l(1)-1) * x2^l(2) * x3^(l(3)-1);
        end
    end
    
elseif axis1 == 2
    if axis2 == 1
        for i=1:N
        l = index(i,:);
        DDPhi3D(i,:) = l(1) * l(2) * x1^(l(1)-1) * x2^(l(2)-1) * x3^l(3);
        end
    elseif axis2 == 2
        for i=1:N
        l = index(i,:);
        DDPhi3D(i,:) = l(2) * (l(2)-1)  * x1^l(1) * x2^(l(2)-2) * x3^l(3);
        end
    elseif axis2 == 3
        for i=1:N
        l = index(i,:);
        DDPhi3D(i,:) = l(2) * l(3) * x1^l(1) * x2^(l(2)-1) * x3^(l(3)-1);
        end
    end 
    
elseif axis1 == 3
    if axis2 == 1
        for i=1:N
        l = index(i,:);
        DDPhi3D(i,:) = l(1) * l(3) * x1^(l(1)-1) * x2^l(2) * x3^(l(3)-1);
        end
    elseif axis2 == 2
        for i=1:N
        l = index(i,:);
        DDPhi3D(i,:) = l(2) * l(3) * x1^l(1) * x2^(l(2)-1) * x3^(l(3)-1);
        end
    elseif axis2 == 3
        for i=1:N
        l = index(i,:);
        DDPhi3D(i,:) = l(3) * (l(3)-1) * x1^l(1) * x2^l(2) * x3^(l(3)-2);
        end
    end 
end
end