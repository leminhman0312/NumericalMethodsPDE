clear all 
clc
clf

jmax = 128;
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
tolCG = 5;

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
u0 = u;

r0 = f-A*u;

%Solving Mz = r 

% z = M\r;  %Do this with CG

%guest z 

z0 = zeros(imax,1);
z0 = cgfunc(M,r0,tolCG,z0,maxiters); %this is z0

p0 = z0;


%% START CHOLESKY

%guest solution u for Au = f

u = zeros(imax,1);

r_cholesky = f-A*u;

err_v = norm(r_cholesky);

for k = 1:maxiters
  
  w = A*p0;
  alpha = ((z0.')*(r_cholesky))/(p0.'*w);
  unew = u + alpha*p0;
  
  rnew_cholesky = r_cholesky-(alpha*w);
  u = unew;
  err_v = [err_v, norm(rnew_cholesky)];
  
  if (norm(rnew_cholesky) < tol)
    break
  end
  
  zold = z0;
  z0 = cgfunc(M,rnew_cholesky,tolCG,zeros(imax,1),maxiters);
  
  
  pnew = z0 + p0*((rnew_cholesky.')*(z0))/(r_cholesky.'*zold);
  
  r_cholesky = rnew_cholesky;
  
  p0 = pnew; 


end  


%% ANALYTICAL 

  
u_analytical = ((xspace.^2)/2)-((xspace.^4)/6)-(1/3);


%for convergence analysis 

umuk = u_analytical-unew.'; %u - uk 

umu0 = u_analytical-u0.';

normuksq = umuk*(A*M)*umuk.';
normu0sq = umu0*(A*M)*umu0.';

lamdamax = max(eig(A*M));
lamdamin = min(eig(A*M));

ratio_lamda = lamdamax/lamdamin;

ratio_norm = normuksq/normu0sq;

rhs_error = 2*((sqrt(ratio_lamda)-1)/(sqrt(ratio_lamda)+1))^(jmax);

lhs_error = ratio_norm;





%% PLOTTING 

figure(1)
plot(xspace,u_analytical,'-r');
title('Incomplete Cholesky Decomposition','FontSize',24)
xlabel('X','FontSize',24)
ylabel('U','FontSize',24)
hold on
plot(xspace,unew,'xk');
legend('ANALYTICAL','NUMERICAL')
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)



figure(2)
plot(log10(err_v).');
title('Cholesky: Error vs iteration number','FontSize',24);
xlabel('Iteration','FontSize',24)
ylabel('Error','FontSize',24)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)


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