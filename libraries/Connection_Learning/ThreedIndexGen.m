function ThreeDIndex = ThreedIndexGen(n)
N = (n+1)*(n+2)/2;
ThreeDIndex = zeros(N,3);
k=1;
for i=1:n+1
    r = i-1;
    Two = TwodIndexGen(n-r);
    [m, l] = size(Two);
    Three = zeros(m,3);
    Three(:,1) = r * ones(m,1);
    Three(:,2:3) = Two;
    ThreeDIndex(k:k+m-1,:)=Three;
    k = k+m;
end
end