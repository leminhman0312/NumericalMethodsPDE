clear all
clc
% clf

%FORM:  uxx = -f = g

%DEFINE GRID SIZE 
a =0; %first grid point 
b = 1; %last grid point 
jmax = 128; %max number of cells
dx = (b-a)/(jmax-1); %spacing 


%DEFINE BASIC PARAMETERS 
tol = 1e-5; %error
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

u_analytical = ((x.^2)/2)-((x.^4)/6)-(1/3)*x;

A = full(gallery('tridiag',jmax,1,-2,1));

A = (1./dx^2)*A;

% 
rhs = (dx^2).*(g);
% 
ureal = (A\(rhs));
ureal = ureal.';

%getting diagonal vectors

D = eye(jmax)*(-2/dx.^2);
Lt = tril(A);
L = Lt - D; 
% 

omega = 0.90;

I = eye(jmax);

u_jacobi = zeros(jmax,1);
u_updated = u_jacobi;
Dinv = inv(D);

res = g - (A*u_jacobi); 

err_v = norm(res);

err_start  = norm(res);


for k = 1:100000
    u_updated = omega*((I-Dinv*A)*u_jacobi+Dinv*g)+(1-omega)*u_jacobi;
    u_jacobi = u_updated; 
    res = g-(A*u_updated);
    err_v = [err_v,norm(res)];
    if (norm(res)/norm(g) < tol)
        break;
    end
    
end


%% CONVERGENCE ANALYSIS

% Binv = Dinv;
% 
% R_jacobi = eye(jmax)-Binv*A;
% max(abs(eig(R_jacobi)))
% R_iterative = omega*(eye(jmax)-Binv*A)+(1-omega)*eye(jmax);
% disp(size(R_iterative));
% Rtrans = R_iterative.';
% norm_R_iter = norm(R_iterative);
% disp(norm_R_iter)
% 
% %spectral radius 
% 
% spec_rad = max(abs(eig(R_iterative)));
% disp(spec_rad);


%% PLOTTING 

f1 = figure(1); 
plot(x,u_analytical,'-r')
hold on
plot(x,u_jacobi,'xk')
title(['Relaxed Jacobi, converge at k = ' num2str(k) ' w = ' num2str(omega)],'FontSize',12)
xlabel('X','FontSize',12)
ylabel('U','FontSize',12)
xt = get(gca, 'XTick');
legend('ANALYTICAL','NUMERICAL')
saveas(f1,'jacobi75.jpg')

% 
f2 = figure(2);
plot((err_v));
title('Error vs iteration number','FontSize',12);
xlabel('Iteration','FontSize',24)
ylabel('Error','FontSize',24)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
saveas(f2,'jacobi_error.jpg')