% Simulating the inviscid Burgers' equation (2-D Convection) by the Finite Difference Method
...(a time march)
% Numerical scheme used is a first order upwind in both space and time

%%
%Specifying Parameters
nx=40;                           %Number of steps in space(x)
ny=40;                           %Number of steps in space(y)       
nt=50;                           %Number of time steps 
dt=0.01;                         %Width of each time step
dx=2/(nx-1);                     %Width of space step(x)
dy=2/(ny-1);                     %Width of space step(y)
x=0:dx:2;                        %Range of x(0,2) and specifying the grid points
y=0:dy:2;                        %Range of y(0,2) and specifying the grid points
u=zeros(nx,ny);                  %Preallocating u
un=zeros(nx,ny);                 %Preallocating un
v=zeros(nx,ny);                  %Preallocating v
vn=zeros(nx,ny);                 %Preallocating vn

%%
% Initial Conditions
for i=1:nx
    for j=1:ny
        if ((0.5<=y(j))&&(y(j)<=1)&&(0.5<=x(i))&&(x(i)<=1))
            u(i,j)=4;
            v(i,j)=1;
        else
            u(i,j)=1;
            v(i,j)=2;
        end
    end
end

%%
%Boundary conditions
u(1,:)=0;
u(nx,:)=0;
u(:,1)=0;
u(:,ny)=0;
v(1,:)=0;
v(nx,:)=0;
v(:,1)=0;
v(:,ny)=0;
i=2:nx-1;
j=2:ny-1;

%%
%Explicit method with F.D in time and B.D in space
for it=0:nt
    un=u;
    vn=v;
    h=quiver(x,y,u',v','Color','black');       %plotting the velocity field
    axis([0 2 0 2])
    title({'2-D Convection';'Transport property vector field {\bfu}=(u_x,u_y)';['time(\itt) = ',num2str(dt*it)]})
    xlabel('Spatial co-ordinate (x) \rightarrow')
    ylabel('Spatial co-ordinate (y) \rightarrow')
    drawnow; 
    refreshdata(h)
    u(i,j)=un(i,j)-(dt*un(i,j).*(un(i,j)-un(i-1,j))/dx)-(dt*vn(i,j).*(un(i,j)-un(i,j-1))/dy);
    v(i,j)=vn(i,j)-(dt*un(i,j).*(vn(i,j)-vn(i-1,j))/dx)-(dt*vn(i,j).*(vn(i,j)-vn(i,j-1))/dy);
    %Boundary Conditions
    u(1,:)=0;
    u(nx,:)=0;
    u(:,1)=0;
    u(:,ny)=0;
    v(1,:)=0;
    v(nx,:)=0;
    v(:,1)=0;
    v(:,ny)=0;
end