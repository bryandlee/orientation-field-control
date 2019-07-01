function DPhi3D = dphi3D(x,axis,index)
x1 = x(1);
x2 = x(2);
x3 = x(3);
[N, b] = size(index);
DPhi3D = zeros(N,1);
if axis == 1
for i=1:N
    l = index(i,:);
    DPhi3D(i,:) = l(1) * x1^(l(1)-1) * x2^(l(2)) * x3^(l(3));
end
elseif axis == 2
for i=1:N
    l = index(i,:);
    DPhi3D(i,:) = l(2) * x1^(l(1)) * x2^(l(2)-1) * x3^(l(3));
end    
elseif axis == 3
for i=1:N
    l = index(i,:);
    DPhi3D(i,:) = l(3) * x1^(l(1)) * x2^(l(2)) * x3^(l(3)-1);
end
end
end