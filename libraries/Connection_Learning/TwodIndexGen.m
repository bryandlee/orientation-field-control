function TwoDIndex = TwodIndexGen(n)
TwoDIndex = zeros(n+1,2);
TwoDIndex(1,:)= [0, n];
for i=1:n
    r = i+1;
    TwoDIndex(r,:)=TwoDIndex(r-1,:) + [1, -1];
end
end