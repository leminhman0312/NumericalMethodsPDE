clear all 
clc

jmax = 5;
imax = jmax;
A = full(gallery('tridiag',jmax,1,-2,1));

A(1,1) = -1; 
A(1,2) = 1;
xspace = linspace(0,1,jmax); %mesh space 
f = 1-2*(xspace'.^2);
dx = 1/(jmax-1);
maxiters = 1000;

%define tolerance 
tol = 1e-4;


%Analytical 

u_analytical = ((xspace.^2)/2)-((xspace.^4)/6)-(1/3);

u = zeros(jmax,1);

f = -f;



%Form Incomplete Cholesky Matrix L 


A = (1/dx^2)*(-A);
L = zeros(imax,jmax);


for i = 1: imax
  sum1 = 0; 
  
  for k = 1:i-1
    sum1 = sum1+ L(i,k)^2;
    
  end
  L(i,i) = (A(i,i) - sum1)^(0.5);
  for j = i+1:jmax
    sum2 = 0;
    for k = 1:i-1
      sum2 = sum2 + (L(i,k).*L(j,k));
      
    end
    
    L(j,i) = (1./(L(i,i)))*(A(j,i) -sum2);
  end
end

    




%calculate M

M = L*L.';
u = zeros(imax,1);

r = f-A*u;

%Solving Mz = r 

% z = M\r;  %Do this with CG


z0 = cgfunc(M,f,1e-4,u,maxiters); %this is z0




p0 = z0;


%% START CHOLESKY

%guest solution u for Au = f

u = zeros(imax,1);

r_cholesky = f-A*u;

err_v = [];

for k = 1:maxiters
  w = A*p0;
  alpha = ((z0.')*(r_cholesky))/(p0.'*w);
  unew = u + alpha*p0;
  
  
  rnew_cholesky = r_cholesky-(alpha*w);
  
  if (norm(rnew_cholesky) < tol)
    r_cholesky = rnew_cholesky;
    u = unew;
    err_v = rnew_cholesky;
    break
    
  zold = z0;  
  znew = cgfunc(M,f,1e-2,u,maxiters);
  
  
  pnew = znew + ((rnew_cholesky.')*(znew))/(r_cholesky.'*zold);
  
  p0 = pnew; 
  
  
  end


end





  


%% ANALYTICAL 

  
u_analytical = ((xspace.^2)/2)-((xspace.^4)/6)-(1/3);





%% PLOTTING 

% plot(xspace,u_analytical,'-r');
% hold on
% plot(xspace,unew,'xk');








function h = cgfunc(A,f,tol,u,maxiters)
  r = f-A*u;
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
      break      
    end      

    beta = (rnewt*rnew)/(rt*r);    
    
    p = rnew+beta*p;
    
    %update value 
    r = rnew;
    u = unew;
%     p = pnew;
  end
  
  
  
  h = unew;
  
  
end  











































%%     CHOLESKY ALGORITHM
%{
for k = 1:1000

  w = A*p;
  alpha = ((z.')*(r))/(p.'*w);
  unew = u + alpha*p;
  rold = r;
  zold = z;
  rnew = r - alpha*w;
  if norm(rnew)<tol
    disp(rnew)
    r = rnew;
    u = unew;
    break
  end
 
  
  znew = M\rnew;
  
  pnew = z + ((rnew.'*znew)/(rold.'*zold))*p;
  
  p = pnew;
  
  
end


plot(xspace,unew);
%}
