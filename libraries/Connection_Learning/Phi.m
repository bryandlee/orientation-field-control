function phi = Phi( x,order )
%PHI 이 함수의 요약 설명 위치
%   자세한 설명 위치
s = size(x);
if s(1) == 2
    n= (order +1)*(order + 2)/2;
    phi = zeros(n,1);
    k=1;
    for i=1:order+1
        for j=1:i
            phi(k,1) = x(1)^(i-j)*x(2)^(j-1);
            k=k+1;
        end
    end
    phi=reshape(phi,[k-1,1]);
elseif s(1) == 3
    n= (order^3 + 6 * order^2 + 11 * order + 6)/6;
    phi = zeros(n,1);
    k=1;
    for i=1:order+1
        for l=1:i
            for j=1:(i-l+1)
                phi(k,1) = x(1)^(l-1) * x(2)^(j-1)  * x(3)^((i-l+1)-j);
                k=k+1;
            end
        end
    end
    phi=reshape(phi,[k-1,1]);
end

end

