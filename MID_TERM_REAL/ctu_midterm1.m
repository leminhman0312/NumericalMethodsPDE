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

% for i = 1:imax+1
%     for j = 1:jmax
%         u(i,j) = 2*y(j);        
%     end
% end
% 
% 
% for i = 1:imax
%     for j = 1:jmax+1
%         
%         v(i,j) = -2*x(i);
%     end
% end

for i = 1: imax
  for j = 1: jmax
    u(i,j) = 2*y(j);
    v(i,j) = -2*x(i);
  end
end







%define initial conditions Q

q = zeros(imax,jmax);

for i = 1:imax
    for j = 1:jmax
        if ((x(i)<0.6)&(x(i)>0.1)&(y(j)>-0.25)&(y(j)<0.25))
            q(i,j) = 1;
        elseif sqrt((x(i)+0.45).^2 + y(j).^2) < 0.35
            q(i,j) = 1-(sqrt((x(i)+0.45).^2 + y(j).^2))...
            /(0.35);
        else
            q(i,j) = 0;
        end
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CORNER TRANSPORT UPWINDING METHOD
%(REMEMBER TO COMMENT OUT DCU)

qnext_ctu = zeros(imax,jmax);



%A and B
aplus = zeros(imax,jmax);
aminus = zeros(imax,jmax);

bplus = zeros(imax,jmax);
bminus = zeros(imax,jmax);


for n = 1: 200
    %define corrective fluxes, zero them out 
    
    disp(n);
   
    g_tilda = zeros(imax,jmax);
    f_tilda = zeros(imax,jmax);
    
    for i = 2:imax-2
        for j = 2:jmax-2         
          

          %DCU fluxes

          up = max(0,u(i,j));
          um = min(0,u(i,j));
          vp = max(0,v(i,j));      
          vm = min(0,v(i,j));

          aplus(i,j) = up*(q(i,j)-q(i-1,j));
          aminus(i,j) = um*(q(i+1,j)-q(i,j));

          bplus(i,j) = vp*(q(i,j)-q(i,j-1));
          bminus(i,j) = vm*(q(i,j+1)-q(i,j));      
          

          qminus_x  = q(i-1,j);
          qminus_x(1,j) = q(imax-1,j);
          
          qplus_x = q(i,j);
          qplus_x(imax,j) = q(1,j);
          
          
          qminus_y = q(i,j-1);
          qminus_y(i,1) = q(i,jmax-1);
          
          qplus_y = q(i,j);
          qplus_y(i,jmax) = q(i,1);
          
          

           %test fluxes
            
          g_tilda(i-1,j) = g_tilda(i-1,j)-0.5*(dt/dx)*(min(0,v(i-1,j)))*(min(0,u(i,j)))*(q(i,j)-q(i-1,j));
          g_tilda(i-1,j+1) = g_tilda(i-1,j+1) -0.5*(dt/dx)*(max(0,v(i-1,j+1)))*(min(0,u(i,j)))*(q(i,j)-q(i-1,j));
          g_tilda(i,j) = g_tilda(i,j)-0.5*(dt/dx)*(min(0,v(i,j)))*(max(0,u(i,j)))*(q(i,j)-q(i-1,j));
          g_tilda(i,j+1) = g_tilda(i,j+1) -0.5*(dt/dx)*(max(0,v(i,j+1)))*(max(0,u(i,j)))*(q(i,j)-q(i-1,j));


          f_tilda(i,j-1) = f_tilda(i,j-1)-0.5*(dt/dy)*(min(0,u(i,j-1)))*(min(0,v(i,j)))*(q(i,j)-q(i,j-1));
          f_tilda(i+1,j-1) = f_tilda(i+1,j-1) -0.5*(dt/dy)*(max(0,u(i+1,j-1)))*(min(0,v(i,j)))*(q(i,j)-q(i,j-1));
          f_tilda(i,j) = f_tilda(i,j)-0.5*(dt/dy)*(min(0,u(i,j)))*max(0,v(i,j))*(q(i,j)-q(i,j-1));
          f_tilda(i+1,j) = f_tilda(i+1,j)-0.5*(dt/dy)*(max(0,u(i+1,j)))*(max(0,v(i,j)))*(q(i,j)-q(i,j-1));        
            

        end
    end
    
    
    %SUM TERMS UP
    
    for i = 2:imax-2
        for j = 2: jmax-2

            %calculate Q 
            
            qnext_ctu(i,j) = q(i,j)-(dt/dx)*(aplus(i,j)+aminus(i,j))-(dt/dy)*(bplus(i,j)+bminus(i,j))...
                -(dt/dx)*(f_tilda(i+1,j)-f_tilda(i,j)) - (dt/dy)*(g_tilda(i,j+1)-g_tilda(i,j));     

                      
        
        end
    end
    
    
    q = qnext_ctu;
    
    figure(1)   
    contourf(x,y,transpose(q))
    colormap(jet)
    colorbar('location','southoutside')
    title(['2D ADVECTION using CTU method at t = ' num2str(n*dt) ' sec'],'Fontsize',20)
    xlabel('X','Fontsize',14,'Fontweight','bold')
    ylabel('Y','Fontsize',14,'Fontweight','bold')
    zlabel('Q','Fontsize',14,'Fontweight','bold')
    
    
    figure(2)    
    surf(x,y,transpose(q))
    view(-37.5,45);
    colormap(jet)
    colorbar('location','southoutside')
    title(['2D ADVECTION using CTU method at t = ' num2str(n*dt) ' sec'],'Fontsize',20)
    xlabel('X','Fontsize',14,'Fontweight','bold')
    ylabel('Y','Fontsize',14,'Fontweight','bold')
    zlabel('Q','Fontsize',14,'Fontweight','bold')

    



end

