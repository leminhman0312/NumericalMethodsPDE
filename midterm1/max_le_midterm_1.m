clear all 
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%define basic properties 
xmin = -1.0;
xmax = 1.0;
ymin = -1.0;
ymax = 1.0;
dx = 1/64;
dy = dx;
imax = ((xmax-xmin)/(dx))+1;
jmax = ((ymax-ymin)/(dy))+1;
tmin = 0.0;
tmax = pi();
dt = 0.4*dx;
nmax = round((tmax-tmin)/dt); %max time step 

%define x and y arrays
x = linspace(xmin,xmax,imax);
y = linspace(ymin,ymax, jmax);

%define u and v matrices
u = linspace(imax+1,jmax);
v = linspace(imax,jmax+1);

for i = 1:imax
    for j = 1:jmax
        u(i,j) = 2*y(j);
        v(i,j) = -2*x(i);
    end
end


%define initial conditions Q

q = zeros(imax,jmax);

for i = 1:imax
    for j = 1:jmax
        q(i,j) = 0.;
        if ((x(i) < 0.6) & (x(i) > 0.1) & (y(j)>-0.25) & (y(j) < 0.25))
            q(i,j) = 1.;
        end
        if sqrt((x(i)+0.45).^2 + y(j).^2) < 0.35
            q(i,j) = 1.-((sqrt((x(i)+0.45)^2 + y(j)^2))...
            /(0.35));
        end
    end
end


% contourf(x,y,transpose(q))



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DONOR CELL UPWINDING METHOD 

qnext = zeros(imax,jmax);
% h = figure;
% filename = '0timestep.gif';
for n = 1:nmax
    disp(n)
    
    qnext = q;
    for i = 2:imax-2
        for j = 2:jmax-2
            
            %calculate u+/- and v+/-
            
%             uplus = 0.5*(u(i,j)+u(i-1,j));
%             uminus = 0.5*(u(i+1,j)+u(i,j));
%             
%             vplus = 0.5*(u(i,j)+u(i,j-1));
%             vminus = 0.5*(u(i,j+1)+u(i,j));
%             
%             uplus  = max(0,uplus);
%             uminus = min(0,uminus);
%             
%             vplus  = max(0,vplus);
%             vminus = min(0,vminus);
            
            uplus = max(0,u(i,j));
            uminus = min(0,u(i,j));
            
            vplus = max(0,v(i,j));
            vminus = min(0,v(i,j));
            
            
            
            
            
            %the finite difference expressions
            forward_diffX_q = q(i+1,j)-q(i,j);
            backward_diffX_q = q(i,j)-q(i-1,j);
            
            forward_diffY_q = q(i,j+1)-q(i,j);
            backward_diffY_q = q(i,j)-q(i,j-1);        
            
            
            %main equation            
            qnext(i,j) = q(i,j)-(dt/dx)*(uplus*backward_diffX_q+uminus*forward_diffX_q)...
                -(dt/dy)*(vplus*backward_diffY_q+vminus*forward_diffY_q);
            
        end
    end
    q=qnext;
    
    %setting BC 
    
    q(imax,:) = q(1,:);
    
    contourf(x,y,transpose(q))
    colormap(jet)
    colorbar('location','southoutside')
    xlabel('x')
    ylabel('y')
    drawnow
    pause(0.01)
    
%   % Capture the plot as an image 
%   frame = getframe(h); 
%   im = frame2im(frame); 
%   [imind,cm] = rgb2ind(im,256); 
%   % Write to the GIF File 
%   if n == 1 
%       imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
%   else 
%       imwrite(imind,cm,filename,'gif','WriteMode','append'); 
%   end 
%     
    
end
 

disp('Done')




