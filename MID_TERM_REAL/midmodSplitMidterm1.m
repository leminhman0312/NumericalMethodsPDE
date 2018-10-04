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
imax = 129%(xmax-xmin)/(dx);
jmax = 129%(ymax-ymin)/(dy);
tmin = 0.0;
tmax = pi();
dt = 0.4*dx;
nmax = round((tmax-tmin)/dt); %max time step 

%define x and y arrays
% x = linspace(xmin,xmax,imax);
% y = linspace(ymin,ymax, jmax);

x = xmin:dx:xmax;
y = ymin:dy:ymax;

%define u and v matrices
u = zeros(imax+1,jmax);
v = zeros(imax,jmax+1);

for i = 1:imax+1
    for j = 1:jmax
        u(i,j) = 2*y(j);
    end
end

for i = 1:imax
    for j = 1:jmax+1
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



qnum = q;
qp1_x = zeros(imax,jmax);
qm1x = zeros(imax,jmax);

qp1y = zeros(imax,jmax);
qm1y = zeros(imax,jmax);
qstar = zeros(imax,jmax);

I = eye(imax);


for n = 1: 200
  
%   disp(n);
  %loop over all i 
  
  qp1_x = circshift(q,-1,1);
  qm1_x = circshift(q,1,1);
  
  
 
  
  %% FOR U 
  
  for i = 1:1:imax
    for j = 1:1:jmax
      %CALCULATE ALL U PLUS AND MINUS
      
      uplus_I_MINUS = max(0.,u(i,j));
      uminus_I_MINUS = min(0.,u(i,j));
      
      
      uplus_I_PLUS = max(0.,u(i+1,j));
      uminus_I_PLUS = min(0.,u(i+1,j));
      
      
      %Compute "A" in leveque pg 103
      
      A_plus_half = uplus_I_PLUS+uminus_I_PLUS;
      A_minus_half = uplus_I_MINUS+uminus_I_MINUS;
      
      
      %Leveque's equation pg 125, HIGHER ORDER
      %FLUXES
      
      F_I_MINUS_HALF = (uminus_I_MINUS*q(i,j)+uplus_I_MINUS*qm1_x(i,j))+0.5*(abs(A_minus_half))*...
      (1.-(dt/dx)*abs(A_minus_half))*(q(i,j)-qm1_x(i,j));
    
      
      F_I_PLUS_HALF = (uminus_I_PLUS*qp1_x(i,j)+uplus_I_PLUS*q(i,j))+0.5*(abs(A_plus_half))*...
      (1.-(dt/dx)*abs(A_plus_half))*(qp1_x(i,j)-q(i,j));
    
    
    
      %LOWER ORDER FLUXES
      
      F_L_MINUS = (uminus_I_MINUS*q(i,j)+uplus_I_MINUS*qm1_x(i,j));
      F_L_PLUS = (uminus_I_PLUS*qp1_x(i,j)+uplus_I_PLUS*q(i,j));
    
    
    
      
      %CALCULATE THETA pg 114
      
      theta_plus_I_PLUS = (q(i,j)-qm1_x(i,j))./(qp1_x(i,j)-q(i,j));
      theta_minus_I_MINUS = (qp1_x(i,j)-q(i,j))./(q(i,j)-qm1_x(i,j));
      
      %shift to get the other THETAS 
      
      %(shift left to get right)
      
      
      if i == imax
          theta_minus_I_PLUS = (q(2,j)-q(1,j))./(q(1,j)-q(imax,j));
      else
        theta_minus_I_PLUS = (qp1_x(i+1,j)-q(i+1,j))./(q(i+1,j)-qm1_x(i+1,j));
      end
      
      %(shift right to get left)
      
      if i == 1
        theta_plus_I_MINUS = (q(imax,j)-q(imax-1,j))./(q(1,j)-q(imax,j));
        
      else
        theta_plus_I_MINUS = (q(i-1,j)-qm1_x(i-1,j))./(qp1_x(i-1,j)-q(i-1,j));
      
      end
    
      
      %Calculate theta plus and theta minus 
      
      theta_MINUS_X = 0.;
      
      if u(i,j) > 0
        theta_MINUS_X = theta_plus_I_MINUS;
      elseif u(i,j)< 0
        theta_MINUS_X = theta_minus_I_MINUS;
      end
      
      
      theta_PLUS_X = 0.;
       if u(i+1,j) > 0
        theta_PLUS_X = theta_plus_I_PLUS;
       elseif u(i+1,j) <0
        theta_PLUS_X = theta_minus_I_PLUS;
      end
      
      
      %CALCULATE PHI
      
      PHI_I_PLUS = max(0,min(1.,theta_PLUS_X));
      PHI_I_MINUS = max(0,min(1.,theta_MINUS_X));
      
      
      
    %Leveque's equation pg 127 (NEED PHI)
    
    
    
%     F_TOTAL_MINUS = (uminus_I_MINUS*q+uplus_I_MINUS*qm1_x)+PHI_I_MINUS*(F_I_MINUS_HALF-(uminus_I_MINUS*q+uplus_I_MINUS*qm1_x));
%     F_TOTAL_PLUS = (uminus_I_PLUS*qp1_x+uplus_I_PLUS*q)+PHI_I_MINUS*(F_I_PLUS_HALF-(uminus_I_PLUS*qp1_x+uplus_I_PLUS*q));
%     
    

     F_TOTAL_PLUS_X = F_L_PLUS + PHI_I_PLUS*(F_I_PLUS_HALF-F_L_PLUS);
     
     F_TOTAL_MINUS_X = F_L_MINUS + PHI_I_MINUS*(F_I_MINUS_HALF-F_L_MINUS);

     
     
     qstar(i,j) = q(i,j) - (dt/dx)*(F_TOTAL_PLUS_X-F_TOTAL_MINUS_X);     
      
     
      
      
      
    end
  end
  
  
  q = qstar;
  
  
  
  %% SHIFT IN V
  
  
  qp1_y = circshift(q,-1,2);
  qm1_y = circshift(q,1,2);
  
  
  
  
  %% FOR V 
    
    for i = 1:imax
      for j = 1: jmax
        
        %CALCULATE VPLUS AND VMINUS
        
        
        vplus_J_MINUS = max(0.,v(i,j));
        vminus_J_MINUS = min(0.,v(i,j));


        vplus_J_PLUS = max(0.,v(i,j+1));
        vminus_J_PLUS = min(0.,v(i,j+1));
        
        
        %Compute "A" in leveque pg 103
      
        A_plus_half = vplus_J_PLUS+vminus_J_PLUS;
        A_minus_half = vplus_J_MINUS+vminus_J_MINUS;
        
        
        
        %Leveque's equation pg 125, HIGHER ORDER
        %FLUXES
      
        F_J_MINUS_HALF = (vminus_J_MINUS*q(i,j)+vplus_J_MINUS*qm1_y(i,j))+0.5*(abs(A_minus_half))*...
        (1.-(dt/dy)*abs(A_minus_half))*(q(i,j)-qm1_y(i,j));


        F_J_PLUS_HALF = (vminus_J_PLUS*qp1_y(i,j)+vplus_J_PLUS*q(i,j))+0.5*(abs(A_plus_half))*...
        (1.-(dt/dy)*abs(A_plus_half))*(qp1_y(i,j)-q(i,j));
      
      
      
        %LOWER ORDER FLUXES
      
        F_L_MINUS = (vminus_J_MINUS*q(i,j)+vplus_J_MINUS*qm1_y(i,j));
        F_L_PLUS = (vminus_J_PLUS*qp1_y(i,j)+vplus_J_PLUS*q(i,j));
        
        
        
    
        %THETA
        
        
            
        theta_plus_J_PLUS = (q(i,j)-qm1_y(i,j))./(qp1_y(i,j)-q(i,j));
        theta_minus_J_MINUS = (qp1_y(i,j)-q(i,j))./(q(i,j)-qm1_y(i,j));
      
        %shift to get the other THETAS 
      
         %(shift left to get right)
      
      
        if j == jmax
            theta_minus_J_PLUS = (q(i,2)-q(i,1))./(q(i,1)-q(i,jmax));
        else
          theta_minus_J_PLUS = (qp1_y(i,j+1)-q(i,j+1))./(q(i,j+1)-qm1_y(i,j+1));
        end

        %(shift right to get left)
      
        if j == 1
          theta_plus_J_MINUS = (q(i,jmax)-q(i,jmax-1))./(q(i,1)-q(i,jmax));

        else
          theta_plus_J_MINUS = (q(i,j-1)-qm1_y(i,j-1))./(qp1_y(i,j-1)-q(i,j-1));

        end
        
        
        
        %THETAMINUS AND THETA PLUS 
        
       %Calculate theta plus and theta minus 
      
        theta_MINUS_Y = 0.;

        if v(i,j) > 0
          theta_MINUS_Y = theta_plus_J_MINUS;
        elseif v(i,j)< 0
          theta_MINUS_Y = theta_minus_J_MINUS;
        end


        theta_PLUS_Y = 0.;
         if v(i,j+1) > 0
          theta_PLUS_Y = theta_plus_J_PLUS;
         elseif v(i,j+1) <0
          theta_PLUS_Y = theta_minus_J_PLUS;
         end
        
        
        %CALCULATE PHI
      
        PHI_J_PLUS = max(0,min(1.,theta_PLUS_Y));
        PHI_J_MINUS = max(0,min(1.,theta_MINUS_Y));
      
         
        
       F_TOTAL_PLUS_Y = F_L_PLUS + PHI_J_PLUS*(F_J_PLUS_HALF-F_L_PLUS);
     
       F_TOTAL_MINUS_Y = F_L_MINUS + PHI_J_MINUS*(F_J_MINUS_HALF-F_L_MINUS);

     
     
      qstar(i,j) = q(i,j) - (dt/dy)*(F_TOTAL_PLUS_Y-F_TOTAL_MINUS_Y);     
      
        
        
        
      end
    end
    
    
    q = qstar;
    
    
    figure(1)   
    contourf(x,y,transpose(q))
    colormap(jet)
    colorbar('location','southoutside')
    title(['2D ADVECTION using Strang Splitting and MinMod limiter at t = ' num2str(n*dt) ' sec'],'Fontsize',20)
    xlabel('X','Fontsize',14,'Fontweight','bold')
    ylabel('Y','Fontsize',14,'Fontweight','bold')
    zlabel('Q','Fontsize',14,'Fontweight','bold')
    
    
    figure(2)    
    surf(x,y,transpose(q))
    view(-37.5,45);
    colormap(jet)
    colorbar('location','southoutside')
    title(['2D ADVECTION using Strang Splitting and MinMod limiter at t = ' num2str(n*dt) ' sec'],'Fontsize',20)
    xlabel('X','Fontsize',14,'Fontweight','bold')
    ylabel('Y','Fontsize',14,'Fontweight','bold')
    zlabel('Q','Fontsize',14,'Fontweight','bold')
  
  
    

   
    
    
  end
  
  
  
  
  
  
  
  




