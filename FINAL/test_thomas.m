clear all 
clc

imax = 4;
d_coeff = 3;
b_coeff = 1;
a_coeff = 1;
diag = d_coeff*ones(imax,1);
below = b_coeff*ones(imax,1);
above = a_coeff*ones(imax,1);
f = [2 4 6 8];
x = zeros(imax,1);
ans = tri(imax,above,below,diag,f);




function thomas = tri(n,a,b,d,c)
    dprime = zeros(n);
    cprime = zeros(n);
    soln = zeros(n,1);
    
    
    dprime(1) = d(1);
    cprime(1) = c(1);
    
    %forward loop 
    for i = 2:n
        dprime(i) = d(i)-(b(i)*a(i-1))/dprime(i-1);
        cprime(i) = c(i)-(b(i)*cprime(i-1))/dprime(i-1);
    end
    
    soln(n) = cprime(n)/dprime(n);
    
    
    %backward loop 
    for i = n-1:-1:2
        soln(i) = (cprime(i)-(a(i)*soln(i+1)))/(dprime(i));
    end

end

