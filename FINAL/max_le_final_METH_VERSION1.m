clear all 
clc
clf

%basic properties 
deltaX = 1./64;
deltaY = deltaX;
dt = 0.1;
tmax = 1;
xmin = 0;
xmax = 1.;
imax = ((xmax-xmin)/deltaX)+1;
jmax = imax;
alpha = 1.;

x = 0:deltaX:1;
y = 0:deltaY:1;


%solution vectors
u = zeros(imax,jmax);
udumb = zeros(imax,jmax);
dx = dt/deltaX^2; %diffusion number in X
dy = dx; %diffusion number in Y 
alpha = 1;


%% THOMAS CONSTANTS X
%define Thomas vectors 

ax = zeros(imax,1); %above
bx = zeros(imax,1); %below
cx = zeros(imax,1); %rhs
diagonalX = zeros(imax,1); %diagonal

for i = 1:imax
    ax(i) = -dx/2;
    bx(i) = -dx/2;
    diagonalX(i) = (1+dx);
end

%% THOMAS CONSTANTS Y
%define Thomas vectors 

ay = zeros(imax,1); %above
by = zeros(imax,1); %below
cy = zeros(imax,1); %rhs
diagonalY = zeros(imax,1); %diagonal


%fill out thomas vectors 

for i = 1:imax
    ay(i) = -dx/2;
    by(i) = -dx/2;
    diagonalY(i) = (1+dx);
end


%% INPUT INITIAL CONDITION 


for i = 1:imax
    for j = 1:jmax
        u(i,j) = initial(x(i),y(j));        
    end
end

% contour(u)
udumb(:,:) = u(:,:);

u_analytical = zeros(imax,jmax);

%% Begin ADI 

for t = 1:tmax
    
    %ANALYTICAL SOLN
    
    for i = 1:imax
        for j = 1:jmax
            u_analytical(i,j) = uexact(t*dt,x(i),y(j));
        end
    end

    
    % X sweep 
    for j = 2:jmax-1
        for i = 2:imax-1
          
            cx(i) = (1-dx)*u(i,j)+(dx/2)*(u(i,j+1)) + (dx/2)*(u(i,j-1)) + (dt/2)*source(t*dt,x(i),y(j));    
        end
        udumb(:,j) = thomas(ax,diagonalX,cx);       
    end
    
    udumb2 = udumb;
    
    %Y sweep  
    for i = 2:imax-1
        for j = 2:jmax-1
           
            cy(j) = (1-dx)*udumb2(i,j)+(dx/2)*(udumb2(i+1,j)) + (dx/2)*(udumb2(i-1,j)) + (dt/2)*source((t+0.5)*dt,x(i),y(j));       
        end
        udumb(i,:) = thomas(ay,diagonalY,cy);      
    end
    
    u = udumb;
       
    err = u - u_analytical.';
end

figure(1)
surf(u)
colormap jet
colorbar
title('Surf plot: Numerical Solution at t = 1. dx = 1/64', 'FontSize',12)
ylabel('Y','FontSize',12)
xlabel('X','FontSize',12)
saveas(figure(1),'numsurf.png')

figure(2)
contour(u,20)
colormap jet
colorbar
title('Contour plot: Numerical Solution at t = 1. dx = 1/64', 'FontSize',12)
ylabel('Y','FontSize',12)
xlabel('X','FontSize',12)
saveas(figure(2),'numcontour.png')


figure(3)
surf(u_analytical)
colormap jet
colorbar
title('Surf plot: Exact Solution at t = 1. dx = 1/64', 'FontSize',12)
ylabel('Y','FontSize',12)
xlabel('X','FontSize',12)
saveas(figure(3),'exactsurf.png')


figure(4)
contour(u_analytical,20)
colormap jet
colorbar
title('Contour plot: Exact Solution at t = 1. dx = 1/64', 'FontSize',12)
ylabel('Y','FontSize',12)
xlabel('X','FontSize',12)
saveas(figure(4),'exactcontour.png')


figure(5)
contour(abs(err),5,'ShowText','on')
colormap jet
colorbar
title('Contour plot: Error distribution at t = 1. dx = 1/64', 'FontSize',12)
ylabel('Y','FontSize',12)
xlabel('X','FontSize',12)
saveas(figure(5),'contour_error.png')

figure(6)
surf(abs(err))
colormap jet
colorbar
title('Surf plot: Error distribution at t = 1. dx = 1/64', 'FontSize',12)
ylabel('Y','FontSize',12)
xlabel('X','FontSize',12)
saveas(figure(6),'surf_error.png')


function exact = uexact(t,x,y)
  exact = exp(-t).*sin(pi.*x).*sin(pi.*y);  
end
  


function f = source(t,x,y)
    f = exp(-t)*sin(pi.*x).*sin(pi.*y)*(2*(pi^2)-1);
end


function f0 = initial(x,y)
    f0 = sin(pi.*x).*sin(pi.*y);
end



function z = thomas(a,dia,f)
   
    n = length(dia);
    z = zeros(n,1);
    p = zeros(n,1);
    q = zeros(n,1);
    
    q(2) = f(2)./dia(2);
    p(2) = a(2)./(dia(2));
    
    for i = 3:n-1
        denom = dia(i)-a(i)*p(i-1);
        p(i) = a(i)./(denom);
        q(i) = (f(i)-a(i)*q(i-1))./denom;
    end
    
    
    z(n) = q(n);
    for i = n-1:-1:1
        z(i) = -p(i)*z(i+1)+q(i);
    end   
    
    
end





