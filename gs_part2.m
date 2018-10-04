clear all
clc

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
% u_analytical = u_analytical.';


A = full(gallery('tridiag',jmax,1,-2,1));

A(1,1) = -2; 
A(1,2) = 2;

% 
rhs = (dx^2).*(g);
% 
ureal = (A\(rhs));
ureal = ureal.';


% z = inv(A);
% ureal = z*rhs;
% ureal = ureal.';





%STARTING GAUSS-SEIDEL LOOP 
% 
%iteration loop 



%getting diagonal vectors

D = ones(jmax,1)*-2;
L = tril(A);
% 
B = L;
Binv = inv(B); 
% 
R = eye(jmax)-Binv*A;


u_gs = zeros(jmax,1);
c = Binv*(dx^2*g);
for k = 1:50000
  
  u_gs = (R*u_gs)+c;
  res = (dx^2)*g-(A*u_gs);
  if (norm(res)<tol)
    break
  end
  
end

disp(norm(res));












% for iter = 1: 10000
%    
%     %Gauss-Seidel loop 
%     for j = 1:jmax
%         if j == jmax
%           unew(j) = 0;
%         elseif j == 1
%           unew(j) = unew(2);
%         else          
%           unew(j) = ((1/2)*(uold(j+1)+unew(j-1)-((dx.^2)*g(j))));
%         end
%     end
%     
%     
%     
%     err(iter) = norm(unew-uold);
%     resi = (dx^2)*g - A*(unew.');
% %     disp(norm(resi));
% %     disp(iter);
%     if err(iter)<tol
%         break 
%     else
%     end
%     uold = unew;
%     iter = iter+1;
%     
%     
%    
%    
%     
%     
% end


% disp(norm(resi));
% disp(unew(1)); 
% disp(unew(jmax));





% res = abs(ureal-u_analytical);
% disp(sum(res)/jmax);

% plot(x,ureal','r')
% hold on
plot(x,u_analytical,'xk')
hold on
plot(x,u_gs,'sb')


% legend('MATLAB','INTEGRAL','NUMERICAL')

























