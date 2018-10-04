clear all 
clc

js = 16;
jb = 512;

jrange = (jb-js);

j_vector = js:1:jb;
err_store = zeros(1,length(j_vector));

for iter = 1: length(err_store)
% jmax = 128

  
  jmax = j_vector(iter);
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

  u = zeros(jmax,1);
  u0 = zeros(jmax,1);

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
  
  error_num = norm(u_analytical-unew.');
  err_store(iter) = error_num;
end




figure(1)
plot(err_store,'-xr');
title('Cholesky: Error vs Gridsize','FontSize',24);
xlabel('Gridsize','FontSize',24)
ylabel('Error','FontSize',24)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)

