function Index3D = Indenx3DGen(n)
N = (n^3 + 6*n^2 + 11*n + 6)/6;
Index3D = zeros(N,3);
k=1;
for i=1:n+1
    r = i-1;
    Index3D(k:k+(r+1)*(r+2)/2-1,:)=ThreedIndexGen(r);
    k=k+(r+1)*(r+2)/2;
end
end