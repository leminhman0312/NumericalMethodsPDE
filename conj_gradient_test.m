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

p = r;





for k = 1:maxiters

    rt = r.';
    pt = p.';
    w = A*p;
    
        
    alpha = (rt*r)/(pt*w);
    unew = u + alpha.*p;
    rnew = r-alpha*w;
    rnewt = rnew.';
    
    
    if (norm(rnew)<tol)
      r = rnew;
      u = unew;
      disp(rnew)
%       p = pnew;
      break
      
    end
  
      

    beta = (rnewt*rnew)/(rt*r);    
    
    p = rnew+beta*p;
    
    %update value 
    r = rnew;
    u = unew;
%     p = pnew;
end

plot(xspace,unew.','xk');
hold on
plot(xspace,u_analytical, '-r');

