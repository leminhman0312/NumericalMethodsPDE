clear all
clc
xmin = 0;
xmax = 1;
nx = 100;
deltax = (xmax-xmin)/(nx-1);
dt = 0.01;
lamda = 1.;
gamma = lamda*dt/deltax;

alpha = 1.0;

dx = alpha*dt/deltax^2;

minus1plus2dx = (1+dx);

xrange = 0:deltax:1;
u = zeros(1,nx);

tboundary = 300;

u = sin(pi*xrange);

%thomas vectors
a = zeros(1,nx);
b = zeros(1,nx);
c = zeros(1,nx);
d = zeros(1,nx);

d(1:end) = minus1plus2dx;
a(2:end-1) = -dx/2;
b(2:end-1) = -dx/2;





% d(1) = -1.0;
% d(end) = -1.0;


%right handed vector 
t0 = 100;
% 
f(1) = -1*tboundary;
f(end) = -1*tboundary;
f(2:end-1) = -1*t0;


f = zeros(1,nx);

%call thomas

for t = 1:10
    u_analytical = (1+(exp((pi^2-1)*(t*dt))-1)/(pi^2-1))*sin(pi*xrange)*exp(-pi^2*t*dt);
    
    
    
    for i = 2:nx-1
       f(i) = ((1-dx)*u(i)+(dx/2)*(u(i-1))+(dx/2)*(u(i+1)))+ dt*(sin(i*deltax*pi)*exp(-t*dt));

        
    end
    u_update = thomas(a,d,f,u);
    u = u_update;
    clf
    plot(u,'-b')
    hold on
    plot(u_analytical,'xr')
    
    pause(0.1)
end

% plot(u)


%% Thomas Algorithm, a = diag, b = upper, c = below 
function [y] = tridiag(iter,diag,upper,below,f)

n = length(f);
v = zeros(n,1);
y = v;
w = diag(1);
y(1) = f(1)/w;

    %Forward substitution
    for i = 2:n
        v(i-1) = below(i-1)/w;
        w = diag(i) - upper(i)*v(i-1);
        y(i) = (f(i)-upper(i)*y(i-1) ) /w;
    end

    %Backward substitution
    for j = n-1:-1:1
       y(j) = y(j) - v(j)*y(j+1);
    end


end



function [z] = thomasold(above,below,dia,f)
    
    n = length(above);
    dprime = zeros(n,1);
    cprime = zeros(n,1);
    dprime(2) = dia(2);
    cprime(2) = f(2);

    %forward
    
    for i = 3:n
       dprime(i) = dia(i) - ((below(i)*above(i-1))/(dprime(i-1)));
       cprime(i) = f(i) - ((below(i)*cprime(i-1))/(dprime(i-1)));
        
    end

    z(n) = cprime(n)/dprime(n);
    
    
    %backward 
    
    for i = n-1:-1:3
        z(i) = (cprime(i)-(above(i)*z(i+1)))/dprime(i);
    end
    
end

function z = thomas(a,dia,f,soln)
   
    n = length(dia);
    z = zeros(n,1);
    p = zeros(n,1);
    q = zeros(n,1);
    
    q(2) = f(2)./dia(2);
    p(2) = -a(2)./(dia(2));
    
    for i = 3:n-1
        denom = dia(i)+a(i)*p(i-1);
        p(i) = -a(i)./(denom);
        q(i) = (f(i)-a(i)*q(i-1))./denom;
    end
    
    
    z(n) = q(n);
%     z(n) = 300.;
%     z(1) = 300.
    for i = n-1:-1:1
        z(i) = p(i)*z(i+1)+q(i);
    end     
    
    
    
    
    
end




