clear all 
clc

jmax = 128;
A = full(gallery('tridiag',jmax,1,-2,1));

A(1,1) = -2; 
A(1,2) = 2;
xspace = linspace(0,1,jmax); %mesh space 
f = 1-2*(xspace'.^2);
dx = 1/(jmax-1);
maxiters = 1000;

%define tolerance 
tol = 1e-4;


%Analytical 

u_analytical = ((xspace.^2)/2)-((xspace.^4)/6)-(1/3);



% define guess parameters 

u = zeros(128,1);

r = ((dx^2)*f)-A*u;

p0 = r;




err_v = norm(r);

for k = 1:maxiters

    rt = r.';
    pt = p0.';
    w = A*p0;
    
        
    alpha = (rt*r)/(pt*w);
    unew = u + alpha.*p0;
    rnew = r-alpha*w;
    rnewt = rnew.';
    u = unew;
    err_v = [err_v,norm(rnew)];
    if (norm(rnew)<tol)
       break      
    end
  
      

    beta = (rnewt*rnew)/(rt*r);    
    
    pnew = rnew+beta*p0;
    
    %update value 
    r = rnew;
    u = unew;
    p0 = pnew;
end

figure(1)
plot(xspace,u_analytical,'-r');
title('Conjugate Gradient without Preconditioner','FontSize',24)
xlabel('X','FontSize',24)
ylabel('U','FontSize',24)
hold on
plot(xspace,unew,'xk');
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
legend('ANALYTICAL','NUMERICAL')


err_v(1) = NaN;
figure(2)
plot((err_v));
title('CG without preconditioner: Error vs iteration number','FontSize',24);
xlabel('Iteration','FontSize',24)
ylabel('Error','FontSize',24)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)


