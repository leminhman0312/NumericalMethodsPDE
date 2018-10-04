clear all 
clc
%define basic propperties 
%solving qt + aqx = 0 where a = 1
% u = q in leveque

imax = 5; %number of points 
dx = 1/imax; 
dt = 0.5*dx;
a = 0; 
b= 1;

x = linspace(a+0.5*dx,b+0.5*dx,imax);
nmax = 1;

u = zeros(1,imax);

%define initial conditions 

for i = 1:imax
    if(x(i)>=0.6 && x(i) <=0.8)
    u(i) = 1.0;   
    elseif(x(i)>=0.2 && x(i)<=0.4)
    u(i) = sin(5.0*pi*(x(i)-0.2));
    else 
    u(i) = 0.0;
    end
end

realu = u; 
numericalu = u;





%Du's method
for n = 1: nmax
    %split U into Up1 and Um1
%     up1 = [numericalu(2:imax),numericalu(1)];
%     um1 = [numericalu(imax),numericalu(1:imax-1)];
    
    um1 = circshift(u,1);
    up1 = circshift(u,[0,-1]);
    
    %compute fluxes
    
    fl = numericalu; %lower flux
    fh = numericalu+0.5*(1-dt/dx)*(up1-numericalu); %higher flux, LW 
    
    %calculate theta
    
    theta = (numericalu-um1)./(up1-numericalu);
    
    %calculate phi, for Min Mod 
    
    phi = max(0,min(theta,1));
    
    
    %flux at j+1/2
    
    flux_right = fl+phi.*(fh-fl);
    
    %do periodic BC, flux at j-1/2 
    flux_left = [flux_right(imax), flux_right(1:imax-1)];
    
    %calculate new u 
    unew = numericalu-dt/dx*(flux_right-flux_left);
    
    %update u 
    numericalu = unew; 


end

























    

%RAM's stuff

% for n = 1:nmax
%     
%     
%     for i = 2:imax-1
%     %calculate u+ and u-
%     
%     uplus = max(0,u(i));
%     uminus = min(0,u(i));
%     
%     
%     %zero out all fluxes
%     F_lower_minus = zeros(1,imax);
%     F_lower_plus = zeros(1,imax);
%     
%     F_higher_plus = zeros(1,imax);
%     F_higher_minus = zeros(1,imax);
%     
%     %FILL IN FLUXES 
%     
%     %LOWER FLUXES
%     
% %     F_lower_minus = uplus*(u(i-1));
% %     F_lower_plus = uminus*(u(i-1));
% 
%     fl = u;
%     fh = u+0.5*(1-dt/dx)*(
%     
%     %HIGHER FLUXES
%     
% %     F_higher_minus = (uplus*u(i-1)) + 0.5*(uplus*dt/dx)*(u(i)-u(i-1));
% %     F_higher_plus = (uminus*u(i+1)) + 0.5*(uminus*dt/dx)*(u(i+1)-u(i));
%     
%     %CALCULATE THETA
%     
%     theta_plus = zeros(1,imax);
%     theta_minus = zeros(1,imax);
%     theta_plus(i) = (u(i)-u(i-1))/(u(i+1)-u(i));
%     theta_minus(i) = (u(i+1)-u(i))/(u(i)-u(i-1));
%     
%     %CALCULATE PHI
%     
%     %PHI MINUS
%     phi_minus = zeros(1,imax);
%     if (theta_minus(i)>=0) && (theta_minus(i)<=1)
%         phi_minus(i) = theta_minus(i);
%     else
%         phi_minus(i) = 1;
%     end
%     
%     
%     %PHI PLUS 
%     
%     phi_plus = zeros(1,imax); 
%     if (theta_plus(i)>=0) && (theta_plus(i)<=1)
%         phi_plus(i) = theta_plus(i);
%     else
%         phi_minus(i) = 1;
%     end
%     
%     
%     %CALCULATE TOTAL FLUXES 
%     
%     Fminus = F_lower_minus+phi_minus*(F_higher_minus-F_lower_minus);
%     Fplus = F_lower_plus+phi_plus*(F_higher_plus-F_lower_plus);
%     
%     %calculate new U 
%     
%     unew = u - (dt/dx)*(Fplus-Fminus);
%     u = unew;
%     
%     end
%     
% 
% end
% 
% plot(x,u);
% 

