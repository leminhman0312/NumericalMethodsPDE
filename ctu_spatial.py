import numpy as np 
import matplotlib.pyplot as plt 



#declare basic variables 
imax = 32
jmax = 32
xmin = 0
xmax = 1.
ymin = 0
ymax = 1.
dx = (xmax-xmin)/(imax)
dy = (ymax-ymin)/(jmax)
dt = 0.9*dx
tmin = 0
tmax = 4
max_time_step = int((tmax-tmin)/(dt))

#create arrays

x = np.linspace(xmin,xmax,imax)
y = np.linspace(ymin,ymax,jmax)



#setting arbitrary u and v 

u = np.zeros((imax+1,jmax))

v = np.zeros((imax,jmax+1))

q = np.zeros((imax,jmax))


aplus_u= np.zeros_like(q)
aminus_u = np.zeros_like(q)

aplus_v= np.zeros_like(q)
aminus_v = np.zeros_like(q)

bplus_u= np.zeros_like(q)
bminus_u = np.zeros_like(q)

bplus_v= np.zeros_like(q)
bminus_v = np.zeros_like(q)


xp = np.linspace(-0.5*dx, 1. + 0.5*dx, imax +1)



#create empty corrective fluxes 

#setting up inital conditions

for i in range(imax):
	for j in range(jmax):
		if (np.sqrt(((x[i]-0.5)**2)+((y[j]-0.5)**2)) < 0.233):
			u[i][j]= 1.-np.sqrt((x[i]-0.5)**2.+(y[j]-0.5)**2.)/(0.233)

		else: 
			u[i][j] = 0.




# start of the loop 

for n in range(int(max_time_step)):	 

	#zero out dummy u
	q_next = np.zeros_like(q)


	#set all corrective flux to zero 

	f_tida_u = np.zeros((imax+1,jmax))
	g_tida_u = np.zeros((imax,jmax+1))

	f_tida_v = np.zeros((imax+1,jmax))
	g_tida_v = np.zeros((imax,jmax+1))

	q = np.zeros((imax, jmax))
	q_next = np.zeros((imax,jmax))






	#DCU ALGORITHM
	for i in range(1,imax-1):
		for j in range(1,jmax-1):


			#calculate avg 

			uavg = 0.5*(u[i][j]+u[i+1][j]) #q horizontal 
			vavg = 0.5*(v[i][j]+v[i][j+1]) #q vertical


			#calculate LW flux , (cell average)

			


			uplus = np.max([0,u[i][j]])
			uminus = np.min([0,u[i][j]])


			uplusp = np.max([0,u[i][j-1]])
			uminusp = np.min([0,u[i][j-1]])

			vplus = np.max([0,v[i][j]])
			vminus = np.min([0,v[i][j]])

			vplusp = np.max([0,v[i-1][j]])
			vminusp = np.min([0,v[i-1][j]])


			#Calculate Gu

			g_tida_u[i-1][j] += -0.5*(dt/dx)*vminusp*uminus*(0.5*(u[i+1][j]-u[i-1][j]))
			g_tida_u[i-1][j+1] += -0.5*(dt/dx)*vplusp*uminus*(0.5*(u[i+1][j]-u[i-1][j]))
			g_tida_u[i][j] += -0.5*(dt/dx)*vminus*uplus*(0.5*(u[i+1][j]-u[i-1][j]))
			g_tida_u[i][j+1] += -0.5*(dt/dx)*vplus*uplus*(0.5*(u[i+1][j]-u[i-1][j]))


			f_tida_u[i][j-1] += -0.5*(dt/dy)*uminusp*vminus*(0.5*((u[i][j]-u[i][j-1])+(u[i+1][j]-u[i+1][j-1])))
			f_tida_u[i+1][j-1] += -0.5*(dt/dy)*uplusp*vminus*(0.5*((u[i][j]-u[i][j-1])+(u[i+1][j]-u[i+1][j-1])))
			f_tida_u[i][j] += -0.5*(dt/dy)*uminus*vplus*(0.5*((u[i][j]-u[i][j-1])+(u[i+1][j]-u[i+1][j-1])))
			f_tida_u[i+1][j] += -0.5*(dt/dy)*uplus*vplus*(0.5*((u[i][j]-u[i][j-1])+(u[i+1][j]-u[i+1][j-1])))



			#A and B for DCU 

			up = np.max([u[i][j],0]) 
			um = np.min([u[i+1][j],0])

			vp = np.max([v[i][j],0])
			vm = np.min([v[i][j+1],0])

			
			aplus_u[i][j] = up*(u[i+1][j]-u[i-1][j])*0.5
			aminus_u[i][j] = um*(u[i+2][j]-u[i][j])*0.5

			bplus_u[i][j] = vp*0.5*((u[i][j]-u[i][j-1])+(u[i+1][j]-u[i+1][j-1])) 
			bminus_u[i][j] = vm*0.5*((u[i][j+1]-u[i][j])+(u[i+1][j+1]-u[i+1][j])) 

			

	#SUM things up 

	#second  to 2nd last
	for i in range(1,imax):
		for j in range(1,jmax):

			uavg = 0.5*(u[i+1][j]+u[i][j])
			uavg_minus = 0.5*(u[i][j]+u[i-1][j])
			
			u_next = uavg - (dt/dx)*(aplus_u[i][j]+aminus_u[i][j]) - \
				(dt/dy)*(bplus_u[i][j]+bminus_u[i][j])-(dt/dx)*(f_tida_u[i+1][j]-f_tida_u[i][j]) \
				-(dt/dy)*(g_tida_u[i][j+1]-g_tida_u[i][j]) #center values from 1 to imax-1

			u_next_minus = uavg_minus - (dt/dx)*(aplus_u[i-1][j]+aminus_u[i-1][j]) - \
				(dt/dy)*(bplus_u[i-1][j]+bminus_u[i-1][j])-(dt/dx)*(f_tida_u[i][j]-f_tida_u[i-1][j]) \
				-(dt/dy)*(g_tida_u[i-1][j+1]-g_tida_u[i-1][j]) #center values from 0 to imax-2

			u[i][j] = 0.5*(u_next+u_next_minus)	


	#set BC

	u[-1,:] = u[-2,:]

	u[0,:] = u[-1,:]
	




	plt.clf()
	plt.contourf(xp,y,u.T,15)
	plt.draw()
	plt.pause(0.001)

	

			











		
