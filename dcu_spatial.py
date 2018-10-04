import numpy as np 
import matplotlib.pyplot as plt 



#declare basic variables 
imax = 80
jmax = 80
xmin = 0
xmax = 1.
ymin = 0
ymax = 1.
dx = (xmax-xmin)/(imax)
dy = (ymax-ymin)/(jmax)
dt = 0.8*dx
tmin = 0
tmax = 4
max_time_step = int((tmax-tmin)/(dt))

#create arrays

x = np.linspace(xmin,xmax,imax)
y = np.linspace(ymin,ymax,jmax)



#setting arbitrary u and v 

u = np.zeros((imax,jmax))

v = np.zeros((imax,jmax))



#setting up inital conditions

for i in range(imax):
	for j in range(jmax):
		if (np.sqrt(((x[i]-0.5)**2)+((y[j]-0.5)**2)) < 0.233):
			u[i][j]= 1-np.sqrt((x[i]-0.5)**2+(y[j]-0.5)**2)/(0.233)

		else: 
			u[i][j] = 0.0


# start of the loop 

for n in range(int(max_time_step)):	 

	


	u_next = np.zeros_like(u)

	for i in range(1,imax-1):
		for j in range(1,jmax-1):


			#calculate LW flux , (cell average)

			uplus = 0.5*(u[i][j]+u[i-1][j])
			uminus = 0.5*(u[i+1][j]+u[i][j])

			vplus = 0.5*(v[i][j]+v[i][j-1])
			vminus = 0.5*(v[i][j+1]+v[i][j])



			uplus = np.max([0,uplus])
			uminus = np.min([0,uminus])

			vplus = np.max([0,vplus])
			vminus = np.min([0,vminus])

			u_next[i][j] = u[i][j] - (dt/dx)*(uplus*(u[i][j]-u[i-1][j])+uminus*(u[i+1][j]-u[i][j])) 
			-(dt/dy)*(vplus*(u[i][j]-u[i][j-1])+vminus*(u[i][j+1]-u[i][j]));

	
	u = u_next.copy()
	u[-1,:] = u[-2,:]

	u[0,:] = u[-1,:]




	plt.clf()
	plt.contourf(x,y,u.T,15)
	plt.draw()
	plt.pause(0.001)



			

			











		
