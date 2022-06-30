'''
EE2703-Assignment-5
Laplace Equation
Syam SriBalaji T
EE20B136
11/03/22
'''

from pylab import*
from scipy.special import*
from numpy import*
from math import*
from sys import*
from mpl_toolkits.mplot3d.axes3d import*
from matplotlib.pyplot import*
import mpl_toolkits.mplot3d.axes3d as p3

Nx = 25
Ny = 25
radius = 8
Niter = 1500
Phi_Matrix = zeros([Ny,Nx])

x1 = linspace(0,Nx,Nx)
y1 = linspace(0,Ny,Ny)
x2 = linspace(-(Ny-1)/2,(Ny-1)/2,Ny)
y2 = linspace(-(Nx-1)/2,(Nx-1)/2,Nx)
Y,X = meshgrid(y2,x2)

#***********************************************************************************************************************************
#Plot of the initial potential configuration
#Figure: 1

index = where((X**2 + Y**2) <= 8*8)
Phi_Matrix[index] = 1
xlabel('x-axis of the plate')
ylabel('y-axis of the plate')
title('Initial Contour plot of the potential')
contourf(x2,y2,Phi_Matrix,cmap = cm.get_cmap("autumn"))
show()

OldPhi_Matrix = Phi_Matrix.copy()

Errors =  list(0 for c in range(0,Niter))
for n in range(0,Niter):
	Phi_Matrix[1:-1,1:-1] = 0.25*(Phi_Matrix[1:-1,0:-2]+ Phi_Matrix[1:-1,2:]+ Phi_Matrix[0:-2,1:-1]+ Phi_Matrix[2:,1:-1])
	Phi_Matrix[1:-1,0] = Phi_Matrix[1:-1,1]
	Phi_Matrix[1:-1,Ny-1] = Phi_Matrix[1:-1,Ny-2]
	Phi_Matrix[0,0:] = Phi_Matrix[1,0:]
	Phi_Matrix[index] = 1
	Errors[n] = (abs(Phi_Matrix-OldPhi_Matrix)).max()
	OldPhi_Matrix = Phi_Matrix.copy()

#***********************************************************************************************************************************
#Plot for Error vs. Number of Iteration
#Figure: 2

NOI = linspace(0,Niter,Niter)
xlabel('no. of Iteration')
ylabel('Error')
title('no. of Iteration vs. Error')
plot(NOI,Errors,'r')
show()

#***********************************************************************************************************************************
#Semilog plot for Error vs. Number of Iteration
#Figure: 3

xlabel('no. of Iteration')
ylabel('Error')
title('Semilog plot for no. of Iteration vs. Error')
semilogy(NOI[::50],Errors[::50],'b')
show()

#***********************************************************************************************************************************
#Loglog plot for Error vs. Number of Iteration
#Figure: 4

xlabel('no. of Iteration')
ylabel('Error')
title('Loglog plot for no. of Iteration vs. Error')
loglog(NOI[::50],Errors[::50],'ro')
show()

#***********************************************************************************************************************************
#Comparision of Fit-1(plot from Vectorization method-Figure: 2) & Fit-2(plot from lstsq method-error minimization)
#Figure: 5

NOI = linspace(0,Niter,Niter)

x_f = zeros((Niter,2))
x_f[:,0] = NOI
x_f[:,1] = 1

log_y_f = zeros((Niter,1))
for i in range(0,Niter):
	log_y_f[i,0] = log(Errors[i])

B1, logA1 = linalg.lstsq(x_f,log_y_f,rcond = None)[0]
A2 = -1*exp(logA1)
B2 = -1*B1
print('The value of A')
print('%e'%A2)
print('The value of B')
print('%e'%B2)

A=exp(logA1)
B=B1

Errors_newf=list( A*exp(B*k) for k in NOI)

xlabel('no. of Iteration')
ylabel('Error')
title('Comparision of Fit-1 & Fit-2')
semilogy(NOI[::50],Errors[::50],'bo',label = 'Fit-1: Values from Vectorization method')
semilogy(NOI[::50],Errors_newf[::50],'ro',label = 'Fit-2: Values from lstsq method')

legend()
show()

#***********************************************************************************************************************************
#Plotting Cumulative error with Stopping condition
#Figure: 6

Cumlt=0
Errors_k = list((float)(-1*A/B)*exp(B*(n+0.5)) for n in NOI)
Errors_cumlt = list(0 for k in NOI)

for i in range(0,Niter):
	Errors_cumlt[i] = Cumlt + Errors_k[i]
	Cumlt += Errors_k[i]
	
xlabel('no. of Iteration')
ylabel('Cumulative error')
title('Stopping condition- Cumulative error vs. no. of Iteration')
semilogy(NOI[::50],Errors_cumlt[::50],'ro')
show()

#***********************************************************************************************************************************
#3-D Surface Plot of Potential
#Figure: 7

from mpl_toolkits.mplot3d import Axes3D
fig1  =  plt.figure(7)
ax = fig1.add_subplot(111, projection="3d")
xlabel('y-axis')
ylabel('x-axis')
ax.set_zlabel('ϕ')
title('3-D surface plot of the potential')
surf  =  ax.plot_surface(Y,X,Phi_Matrix.T, rstride = 1, cstride = 1,cmap = cm.jet,linewidth = 0,antialiased = False)
show()

#***********************************************************************************************************************************
#Contour Plot of the Potential
#Figure: 8

fig,ax = subplots(1,1)
cp  =  ax.contourf(Y,-X, Phi_Matrix)
plot(index[0]-(Ny-1)/2,index[1]-(Nx-1)/2,'ro')

fig.colorbar(cp)
xlabel('x-axis')
ylabel('y-axis')
title('Final Contour plot of the potential')
show()

#***********************************************************************************************************************************
#Vector plot of the current flow
#Figure: 9

Jx = zeros([Ny,Nx])
Jy = zeros([Ny,Nx])
for i in range(1,Nx-1):
	for j in range(1,Ny-1):
		Jy[i,j] = 0.5*(Phi_Matrix[i,j-1]-Phi_Matrix[i,j+1])
		Jx[i,j] = 0.5*(Phi_Matrix[i+1,j]-Phi_Matrix[i-1,j])

quiver(y1,x1,Jy[::-1,:],Jx[::-1,:])
xlabel('x-axis of the plate')
ylabel('y-axis of the plate')
title('Vector plot of current flow')
plot(index[1],index[0],'ro')
show()

#Thank you!

























