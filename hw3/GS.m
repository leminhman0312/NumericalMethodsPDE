clear all
clc

%FORM:  uxx = -f = g

%DEFINE GRID SIZE 
a =0; %first grid point 
b = 1; %last grid point 
jmax = 128; %max number of cells
dx = (b-a)/(jmax); %spacing 


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

rhs = (dx^2).*(g);

ureal = A\(rhs);





%STARTING GAUSS-SEIDEL LOOP 
% 
%iteration loop 
for iter = 1: 10000
   
    %Gauss-Seidel loop 
    for j = 1:jmax
        if j == jmax
          unew(j) = 0;
        elseif j == 1
          unew(j) = unew(jmax);
        else          
          unew(j) = (1/2)*(uold(j+1)+uold(j-1)-((dx.^2)*g(j)));
        end
    end
    
    
    
    err(iter) = norm(unew-uold,2);
    if err(iter)<tol
        break 
    else
    end
    uold = unew;
    iter = iter+1;
    
    
   
   
    
    
end

plot(x,ureal','-r')
hold on
plot(x,unew,'xb')
% plot(x,u_analytical,'ok')
% legend('Exact using MATLAB','Numerical','Exact INTERGAL','location','southeast');































