clear all
clc
clf

%FORM:  uxx = -f = g

%DEFINE GRID SIZE 
a =0; %first grid point 
b = 1; %last grid point 
jmax = 128; %max number of cells
dx = (b-a)/(jmax-1); %spacing 


%DEFINE BASIC PARAMETERS 
tol = 1e-20; %error
x = linspace(a,b,jmax); %mesh space 
unew = zeros(1,jmax); %solution 
%right hand side vector 
g = 1-2.*(x'.^2);

%initial conditions
u = 0*(x');
u(1) = 0; 
u(jmax) = 0; 

uold = u;

%ANALYTICAL SOLN 

u_analytical = ((x.^2)/2)-((x.^4)/6)-(1/3);

A = full(gallery('tridiag',jmax,1,-2,1));


A(1,1) = -2; 
A(1,2) = 2;

A = (1./dx^2)*A;

% 
rhs = (dx^2).*(g);
% 
ureal = (A\(rhs));
ureal = ureal.';

%getting diagonal vectors

D = ones(jmax,1)*-2;
Lt = tril(A);
% 
B = Lt;
Binv = inv(B); 
% 
R = eye(jmax)-Binv*A;


u_gs = zeros(jmax,1);
c = Binv*(g);
res = g - (A*u_gs);
err_v = norm(res);
for k = 1:50000
  
  u_gs = (R*u_gs)+c;
  res = g-(A*u_gs);
  err_v = [err_v, norm(res)];
  if (norm(res)<tol)
     break
  end
  
end

% disp(norm(res));


%% PLOTTING 

figure(1) 
plot(x,u_analytical,'-r')
hold on
plot(x,u_gs,'xk')
title('Gauss Seidel','FontSize',24)
xlabel('X','FontSize',24)
ylabel('U','FontSize',24)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
legend('ANALYTICAL','NUMERICAL')


figure(2)
plot(log10(err_v));
title('Gauss Seidel: Error vs iteration number','FontSize',24);
xlabel('Iteration','FontSize',24)
ylabel('Error','FontSize',24)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)

























