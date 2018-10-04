clear all 
clc

%basic properties 
deltaX = 1./64;
deltaY = deltaX;
dt = 0.001;
tmax = 1.;
xmin = 0;
xmax = 1.;
imax = (xmax-xmin)/deltaX;
jmax = imax;
alpha = 1.;

x = 0:deltaX:1;
y = 0:deltaY:1;


%solution vectors
uhalf = zeros(imax,jmax);
u0 = zeros(imax,jmax);
ufinal = zeros(imax,jmax);
u = zeros(imax,jmax);

dx = dt/deltaX^2; %diffusion number in X
dy = dx; %diffusion number in Y 




%% THOMAS CONSTANTS X
%define Thomas vectors 

ax = zeros(imax+1,1); %above
bx = zeros(imax+1,1); %below
cx = zeros(imax+1,1); %rhs
diagonalX = zeros(imax+1,1); %diagonal



%thomas Constants
d1 = 0.5*dx;
d2 = 0.5*dy;
plus2d2 = 1+(2*d2);
plus2d1 = 1+(2*d1);
minus2d2 = 1-(2*d2);
minus2d1 = 1-(2*d1);


%fill out thomas vectors 

for i = 1:imax
    ax(i) = -d1;
    bx(i) = -d1;
    diagonalX(i) = plus2d1;
end

ax(1) = 0.0;
ax(imax) = 0.0;
bx(1) = 0.0;
bx(imax) = 0.0;
diagonalX(1) = 1.0;
diagonalX(imax) = 1.0;


%% THOMAS CONSTANTS Y
%define Thomas vectors 

ay = zeros(imax+1,1); %above
by = zeros(imax+1,1); %below
cy = zeros(imax+1,1); %rhs
diagonalY = zeros(imax+1,1); %diagonal


%fill out thomas vectors 

for i = 1:imax
    ay(i) = -d2;
    by(i) = -d2;
    diagonalY(i) = plus2d2;
end

ay(1) = 0.0;
ay(imax) = 0.0;
by(1) = 0.0;
by(imax) = 0.0;
diagonalY(1) = 1.0;
diagonalY(imax) = 1.0;


%% RHS vector 
cx = 0;
cy = 0;

t0 = 0;

%% INPUT INITIAL CONDITION 


for i = 1:imax
    for j = 1:jmax
        u(i,j) = initial(x(i),y(j));        
    end
end



%% Begin ADI 

for t = 1:tmax
    % X sweep 
    for j = 2:jmax-1
        for i = 2:imax-1
            if (i ==2) 
                cx(i) = d2*u(i,j+1) + (minus2d2*u(i,j)) + (d2*u(i,j-1)) + (d1*t0);
            elseif (i == imax-1)
                cx(i) = d2*u(i,j+1) + (minus2d2*u(i,j)) + (d2*u(i,j-1)) + (d1*t0);
            else
                cx(i) = d2*u(i,j+1) + (minus2d2*u(i,j)) + (d2*u(i,j-1));
            end
            cx(i) = d2*u(i,j+1) + (minus2d2*u(i,j)) + (d2*u(i,j-1)) + (d1*t0);    
        end
        udumb = tridiag(j,diagonalX,ax,bx,cx);       
    end
    %Y sweep  
    for i = 2:imax-1
        for j = 2:jmax-1
            if (j ==2) 
                cy(j) = d1*u(i+1,j) + (minus2d1*u(i,j)) + (d1*u(i-1,j)) + (d2*t0);
            elseif (j == jmax-1)
                cy(j) = d1*u(i+1,j) + (minus2d1*u(i,j)) + (d1*u(i-1,j)) + (d2*t0);
            else
                cy(j) = d1*u(i+1,j) + (minus2d1*u(i,j)) + (d1*u(i-1,j));
            end
%             cy(j) = d1*udumb(i+1,j) + (minus2d1*udumb(i,j)) + (d1*udumb(i-1,j));    
        end
        ufinal = tridiag(i,diagonalY,ay,by,cy);       
    end
    
    
       
    

    t = t + dt;
    
end




%% TEST THOMAS 

%{


n = 4;
a_coeff = 3;
b_coeff = 1;
c_coeff = 1;
a = a_coeff*ones(n,1);
b = b_coeff*ones(n,1);
c = c_coeff*ones(n,1);
rh = [2 4 6 8];
t = tridiag(a,b,c,rh);
A = diag(a,0) + diag(b_coeff*ones(n-1,1),-1) + diag(c_coeff*ones(n-1,1),1); 

%}











%% Thomas Algorithm, a = diag, b = upper, c = below 
function [y] = tridiag(iter,a,b,c,f)

n = length(f);
v = zeros(n,1);
y = v;
w = a(1);
y(1) = f(1)/w;

    %Forward substitution
    for i = 2:n
        v(i-1) = c(i-1)/w;
        w = a(i) - b(i)*v(i-1);
        y(i) = (f(i)-b(i)*y(i-1) ) /w;
    end

    %Backward substitution
    for j = n-1:-1:1
       y(j) = y(j) - v(j)*y(j+1);
    end


end





function exact = uexact(t,x,y)
  exact = exp(-t)*sin(pi.*x)*sin(pi.*y);  
end
  


function rhs = f(t,x,y)
    rhs = exp(-t)*sin(pi.*x)*sin(pi.*y)*(2*pi^2-1);
end


function f0 = initial(x,y)
    f0 = sin(pi.*x)*sin(pi.*y);
end
  
  






