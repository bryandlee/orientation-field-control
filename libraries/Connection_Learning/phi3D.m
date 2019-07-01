function Phi3D = phi3D(x,index)
x1 = x(1);
x2 = x(2);
x3 = x(3);
[N, b] = size(index);
Phi3D = zeros(N,1);
for i=1:N
    l = index(i,:);
    Phi3D(i,:) = x1^(l(1)) * x2^(l(2)) * x3^(l(3));
end
end