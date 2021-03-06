
clear all
clc

% jmax = 512;

j_length = 10;
jmaxrange = zeros(1,j_length); 
for i = 1: j_length
    jmaxrange(i) = 16*i;
end

k_range = zeros(1,j_length);

% jmaxrange = [16,32,64,128,256];
% k_range = zeros(1,5);
for j = 1:length(jmaxrange)
    jmax = jmaxrange(j);

    x = linspace(0,1,jmax);
    x = x.';
    u = zeros(1,jmax);
    u = u.';
    g = 1-2.*(x'.^2);
    g = g.';

    dx = (1-0)/(jmax-1);

    A = full(gallery('tridiag',jmax,1,-2,1));


    A = (1./dx^2)*A;

    u = zeros(jmax,1);
    res = g-(A*u);

    err_v = norm(res);

    tol = 1e-5;

    kmultigrid = 5000;

    err_start = norm(res);
    for k = 1:kmultigrid
        u = vfunct(u,g,8*dx,1,0);
        res = g-(A*u);
    %     disp(norm(res));
        err_v = [err_v, norm(res)];
        if (norm(res)<tol)        
          break;

        end
    end
    k_range(j) = k;
    
end


fitcurve = jmaxrange.*log10(jmaxrange);

fitcurve(1:3) = [];

jmaxrange(1:3) = [];
k_range(1:3) = [];




plot(jmaxrange,k_range,'-r');
hold on
plot(jmaxrange,fitcurve,'-b')
title('Grid Size vs Iteration Number','FontSize',12);
xlabel('Grid Size','FontSize',24)
ylabel('Iteration Number','FontSize',24)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
legend('Numerical','NlogN Fit');













% disp(k);
% 
% err_end = norm(res);
% 
% rho = err_end/err_start;
% 
% disp(rho);
% 
% 
% u_analytical = ((x.^2)/2)-((x.^4)/6)-(1/3)*x;

  

% f1 = figure(1); 
% plot(x,u_analytical,'-r')
% hold on
% plot(x,u,'xk')
% title(['Multigrid V-cycle, k = ' num2str(k)],'FontSize',12)
% xlabel('X','FontSize',24)
% ylabel('U','FontSize',24)
% xt = get(gca, 'XTick');
% set(gca, 'FontSize', 16)
% legend('ANALYTICAL','NUMERICAL')
% saveas(f1,'vcycle.jpg')
% 
% %   
% f2 = figure(2);
% plot((err_v));
% title('Error vs iteration number','FontSize',12);
% xlabel('Iteration','FontSize',24)
% ylabel('Error','FontSize',24)
% xt = get(gca, 'XTick');
% set(gca, 'FontSize', 16)
% saveas(f2,'vcycle_error.jpg')





%% FUNCTION FOR V-CYCLE 

function vcycle = vfunct(u,f,maxdx,xmax,xmin)


s = size(u);
jmax = s(1);
dx = (xmax-xmin)/(jmax-1);
x = linspace(xmin,xmax,jmax);


%% DEFINE BASIC PARAMETERS 
tol = 1e-5; %error

A = full(gallery('tridiag',jmax,1,-2,1));


A = (1./dx^2)*A;



%% RELAX, DO JACOBI 


D = eye(jmax)*(-2/dx.^2);


omega = 0.75;
B =D;
Binv = inv(B); 


% 
R = eye(jmax)-Binv*A;


u_relax = u;
c = Binv*(f);
res = f - (A*u_relax);
err_v = norm(res);


for k = 1:10
  
  u_new = (R*u_relax)+c;
  u_relax = omega*u_new+(1.-omega)*u_relax;
  
  res = f-(A*u_relax); 
end

%% CHECKING 

if (dx >= maxdx)
  vcycle = u_relax;
else
  r2h = zeros(jmax/2,1);
  for i = 1:jmax/2
    r2h(i) = res(2*i);
  end
  u2h = zeros(jmax/2,1);
  
  u2h = vfunct(u2h,r2h,maxdx,1,0);
  
  x2 = linspace(xmin,xmax,jmax/2);
   
  %INTERPOLATION
  interp_result = interp1(x2.',u2h,x.');
  
  u_relax = u_relax + interp_result;
%   disp(size(u_relax));
  
  %% RELAX AGAIN

  res = f - (A*u_relax);
  err_v = norm(res);

  for k = 1:10

    u_new = (R*u_relax)+c;
    u_relax = omega*u_new+(1.-omega)*u_relax;
    res = f-(A*u_relax);

   end   
    
end


vcycle = u_relax;



end