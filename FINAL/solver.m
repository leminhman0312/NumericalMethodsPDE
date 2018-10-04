function [y] = tridiag(a,b,c,f)




n = length(f);
v = zeros(n,1);
y = v;
w = a(1);
y(1) = f(1)/w;


for i = 2:n
    v(i-1) = c(i-1)/w;
    w = a(i) - b(i)*v(i-1);
    y(i) = (f(i)-b(i)*y(i-1) ) /w;
end


for j = n-1:-1:1
   y(j) = y(j) - v(j)*y(j+1);
end


end
