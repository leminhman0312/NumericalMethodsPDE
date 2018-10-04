import numpy as np 
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D


#declare basic variables 
imax = 80
jmax = 80
xmin = 0
xmax = 1.
ymin = 0
ymax = 1.
dx = (xmax-xmin)/(imax)
dy = (ymax-ymin)/(jmax)
dt = 0.5*dx
tmin = 0
tmax = 4
max_time_step = int((tmax-tmin)/(dt))

#create arrays

x = np.linspace(xmin,xmax,imax)
y = np.linspace(ymin,ymax,jmax)
q = np.zeros((imax,jmax))


#setting up inital conditions

for i in range(imax):
	for j in range(jmax):
		if (np.sqrt(((x[i]-0.5)**2)+((y[j]-0.5)**2)) < 0.233):
			q[i][j]= 1-np.sqrt((x[i]-0.5)**2+(y[j]-0.5)**2)/(0.233)

		else: 
			q[i][j] = 0
	




#setting arbitrary u and v 

u = 1.

v = 1.


# start of the loop 

for n in range(int(max_time_step)):	




	
	#loop to calculate 

	uplus = np.max([0,u])
	uminus = np.min([0,u])

	vplus = np.max([0,v])
	vminus = np.min([0,v])


	q_next = np.zeros_like(q)

	for i in range(1,imax-1):
		for j in range(1,jmax-1):

			
			q_next[i][j] = q[i][j] - (dt/dx)*(uplus*(q[i][j]-q[i-1][j])+uminus*(q[i+1][j]-q[i][j])) -(dt/dy)*(vplus*(q[i][j]-q[i][j-1])+vminus*(q[i][j+1]-q[i][j]));

	
	q = q_next.copy()


	
	plt.clf()
	plt.contourf(x,y,q)
	plt.draw()
	plt.pause(0.001)



			

			











		
