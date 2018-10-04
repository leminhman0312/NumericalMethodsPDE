%CREATED BY: MAX LE
%MIDTERM 1 MTH5315 NUMERICAL METHODS FOR PDE
%DATE: 3/15/2018
tic

disp('START DCU METHOD')
%% DONOR CELL UPWINDING METHOD
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

%{
% %PLOTTING INTIAL CONDITION (UNCOMMENT THIS)
figure(1)
contourf(x,y,transpose(q))
xt = get(gca, 'XTick');
set(gca, 'FontSize', 12)
yt = get(gca, 'YTick');
set(gca, 'FontSize', 12)
colormap(jet)
colorbar('location','southoutside')
title(['2D ADVECTION initial conditions, t = 0 sec '],'Fontsize',14)
xlabel('X','Fontsize',14,'Fontweight','bold')
ylabel('Y','Fontsize',14,'Fontweight','bold')
zlabel('Q','Fontsize',14,'Fontweight','bold')

 
    
figure(2)    
surf(x,y,transpose(q))
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
yt = get(gca, 'YTick');
set(gca, 'FontSize', 16)
view(-37.5,45);
colormap(jet)
colorbar('location','southoutside')
title(['2D ADVECTION initial conditions, t = 0 sec '],'Fontsize',14)
xlabel('X','Fontsize',14,'Fontweight','bold')
ylabel('Y','Fontsize',14,'Fontweight','bold')
zlabel('Q','Fontsize',14,'Fontweight','bold')
%}


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
    
    figure(1)
    contourf(x,y,transpose(q))
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 12)
    yt = get(gca, 'YTick');
    set(gca, 'FontSize', 12)
    colormap(jet)
    colorbar('location','southoutside')
    title(['2D ADVECTION using DCU method at t = ' num2str(n*dt) ' sec'],'Fontsize',14)
    xlabel('X','Fontsize',14,'Fontweight','bold')
    ylabel('Y','Fontsize',14,'Fontweight','bold')
    zlabel('Q','Fontsize',14,'Fontweight','bold')
%     pause(0.01)

     
    
    
    
    figure(2)    
    surf(x,y,transpose(q))
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 16)
    yt = get(gca, 'YTick');
    set(gca, 'FontSize', 16)
    view(-37.5,45);
    colormap(jet)
    colorbar('location','southoutside')
    title(['2D ADVECTION using DCU method at t = ' num2str(n*dt) ' sec'],'Fontsize',20)
    xlabel('X','Fontsize',14,'Fontweight','bold')
    ylabel('Y','Fontsize',14,'Fontweight','bold')
    zlabel('Q','Fontsize',14,'Fontweight','bold')


    
end

disp('END DCU METHOD')




%% CONOR TRANSPORT UPWINDING SCHEME

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  START OF CTU CODE %%%%%%%%%%%%%%%%%%%%%%%%%

disp('STARTING CTU')
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


for n = 1: nmax
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
    
    figure(3)   
    contourf(x,y,transpose(q))
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 16)
    yt = get(gca, 'YTick');
    set(gca, 'FontSize', 16)
    colormap(jet)
    colorbar('location','southoutside')
    title(['2D ADVECTION using CTU method at t = ' num2str(n*dt) ' sec'],'Fontsize',20)
    xlabel('X','Fontsize',14,'Fontweight','bold')
    ylabel('Y','Fontsize',14,'Fontweight','bold')
    zlabel('Q','Fontsize',14,'Fontweight','bold')
    
    
    figure(4)    
    surf(x,y,transpose(q))
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 16)
    yt = get(gca, 'YTick');
    set(gca, 'FontSize', 16)
    view(-37.5,45);
    colormap(jet)
    colorbar('location','southoutside')
    title(['2D ADVECTION using CTU method at t = ' num2str(n*dt) ' sec'],'Fontsize',20)
    xlabel('X','Fontsize',14,'Fontweight','bold')
    ylabel('Y','Fontsize',14,'Fontweight','bold')
    zlabel('Q','Fontsize',14,'Fontweight','bold')
    



end

disp('END CTU')







%% STRANG SPLITTING AND MINMOD LIMITER

disp('STARTING MIN MOD')

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
x = linspace(xmin,xmax,imax);
y = linspace(ymin,ymax, jmax);
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


for n = 1: nmax
  
  disp(n);
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
      (1.-(dt/2/dx)*abs(A_minus_half))*(q(i,j)-qm1_x(i,j));
    
      
      F_I_PLUS_HALF = (uminus_I_PLUS*qp1_x(i,j)+uplus_I_PLUS*q(i,j))+0.5*(abs(A_plus_half))*...
      (1.-(dt/2/dx)*abs(A_plus_half))*(qp1_x(i,j)-q(i,j));
    
    
    
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
      
          

     F_TOTAL_PLUS_X = F_L_PLUS + PHI_I_PLUS*(F_I_PLUS_HALF-F_L_PLUS);
     
     F_TOTAL_MINUS_X = F_L_MINUS + PHI_I_MINUS*(F_I_MINUS_HALF-F_L_MINUS);

     
     
     qstar(i,j) = q(i,j) - (dt/2/dx)*(F_TOTAL_PLUS_X-F_TOTAL_MINUS_X);     
      
     
      
      
      
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
      (1.-(dt/2/dx)*abs(A_minus_half))*(q(i,j)-qm1_x(i,j));
    
      
      F_I_PLUS_HALF = (uminus_I_PLUS*qp1_x(i,j)+uplus_I_PLUS*q(i,j))+0.5*(abs(A_plus_half))*...
      (1.-(dt/2/dx)*abs(A_plus_half))*(qp1_x(i,j)-q(i,j));
    
    
    
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
 
    

     F_TOTAL_PLUS_X = F_L_PLUS + PHI_I_PLUS*(F_I_PLUS_HALF-F_L_PLUS);
     
     F_TOTAL_MINUS_X = F_L_MINUS + PHI_I_MINUS*(F_I_MINUS_HALF-F_L_MINUS);

     
     
     qstar(i,j) = q(i,j) - (dt/2/dx)*(F_TOTAL_PLUS_X-F_TOTAL_MINUS_X);     
      
     
      
      
      
    end
  end   
    
    
    q = qstar;
    
    
    figure(5)   
    contourf(x,y,transpose(q))
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 16)
    yt = get(gca, 'YTick');
    set(gca, 'FontSize', 16)
    colormap(jet)
    colorbar('location','southoutside')
    title(['2D ADVECTION using Strang Splitting and MinMod limiter at t = ' num2str(n*dt) ' sec'],'Fontsize',20)
    xlabel('X','Fontsize',14,'Fontweight','bold')
    ylabel('Y','Fontsize',14,'Fontweight','bold')
    zlabel('Q','Fontsize',14,'Fontweight','bold')
   
    
    figure(6)    
    surf(x,y,transpose(q))
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 16)
    yt = get(gca, 'YTick');
    set(gca, 'FontSize', 16)
    view(-37.5,45);
    colormap(jet)
    colorbar('location','southoutside')
    title(['2D ADVECTION using Strang Splitting and MinMod limiter at t = ' num2str(n*dt) ' sec'],'Fontsize',20)
    xlabel('X','Fontsize',14,'Fontweight','bold')
    ylabel('Y','Fontsize',14,'Fontweight','bold')
    zlabel('Q','Fontsize',14,'Fontweight','bold')
  


   
    
    
end
  disp('ENDING MIN MOD')
   
  toc
  
  
 
  
  
  
  
















