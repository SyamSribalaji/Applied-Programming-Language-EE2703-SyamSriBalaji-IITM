'''
EE2703_Jan-May,2022
End Semester Exam
Syam SriBalaji T
EE20B136
12/05/22
'''

from pylab import*

l=0.5                                                          #Quater wavelength in metre
c=2.9979e8                                                     #Speed of light
mu0=4e-7*pi                                                    #Permeability of free space
N=4                                                            #NO. of section in each half of antenna
Im=1.0                                                         #Current in antenna in ampere
a=0.01                                                         #Radius in metre
lamda=l*4.0                                                    #Wavelength
k=2*pi/lamda      	                                       #Frequency
dz=l/N                                                         #Spacing of current samples

#**********************************************************************************************************
#Question-1

I= zeros(2*N+1)                                                #Calculating current through Standard analysis
z = linspace(-N*dz,N*dz,N*2+1)

I[N:]= sin(k*(l-z[N:]))                                        #Usage of given formula
I[0:N]=I[N+1:][::-1]

#**********************************************************************************************************
#Question-2

p=array(range(0,2*N-2))                                        #Creating a Whole number [0,2N-3] to use instead of for loops 
p=p.reshape((len(p),1))

J = zeros([2*N-2,1])                                           #Defining M & J matrices
M=(1/(2*pi*a))*identity(2*N-2)

#**********************************************************************************************************
#Question-3

P_B = zeros([2*N-2,1],dtype=complex)                           #Defining P_B and P_iN
R_iN = zeros([2*N-2,1])
                                                               #Finding R_iN
R_iN[0:N-1]=((N-1-p[0:N-1])*dz)**2                             #Dividing the steps in formula to prevent error and improve speed
R_iN[0:N-1]+=a**2
R_iN[0:N-1]=sqrt(R_iN[0:N-1])

R_iN[N-1:2*N-2]=((p[N-1:2*N-2]-N+2)*dz)**2
R_iN[N-1:2*N-2]+=a**2
R_iN[N-1:2*N-2]=sqrt(R_iN[N-1:2*N-2])

                                                               #Finding P_B from the formula without for loop
P_B[0:2*N-2]=(mu0/(4*pi))*((exp(-1j*k*R_iN[0:2*N-2]))/R_iN[0:2*N-2])*dz

                                                               #Defining and finding P_ij & R_ij with for loop		
P_ij = zeros([2*N-2,2*N-2],dtype=complex)
R_ij = zeros([2*N-2,2*N-2])

for i in range(0,2*N-2):
	for j in range(0,2*N-2):
		if i <= N-2:
			if j <= N-2:
				R_ij[i][j]=sqrt(a**2+(abs(i-j)*dz)**2)
				P_ij[i][j]=(mu0/(4*pi))*((exp(-1j*k*R_ij[i][j]))/(R_ij[i][j]))*dz
			elif j >= N-1:
				R_ij[i][j]=sqrt(a**2+((abs(i-j)+1)*dz)**2)
				P_ij[i][j]=(mu0/(4*pi))*((exp(-1j*k*R_ij[i][j]))/(R_ij[i][j]))*dz
		elif i >= N-1:
			if j <= N-2:
				R_ij[i][j]=sqrt(a**2+((abs(i-j)+1)*dz)**2)
				P_ij[i][j]=(mu0/(4*pi))*((exp(-1j*k*R_ij[i][j]))/(R_ij[i][j]))*dz
			elif j >= N-1:
				R_ij[i][j]=sqrt(a**2+(abs(i-j)*dz)**2)
				P_ij[i][j]=(mu0/(4*pi))*((exp(-1j*k*R_ij[i][j]))/(R_ij[i][j]))*dz


                                                               #Tried to find R_ij without for loop							
																										
'''
R_ij[0:N-1][0:N-1]=sqrt(a**2+(abs(X[0:N-1][0:N-1]-Y[0:N-1][0:N-1])*dz)**2)
R_ij[0:N-1][N-1:2*N-2]=sqrt(a**2+((X[N-1:2*N-2][0:N-1]-Y[N-1:2*N-2][0:N-1])*dz)**2)
R_ij[N-1:2*N-2][0:N-1]=sqrt(a**2+((X[0:N-1][N-1:2*N-2]-Y[0:N-1][N-1:2*N-2])*dz)**2)
R_ij[N-1:2*N-2][N-1:2*N-2]=sqrt(a**2+(abs(X[N-1:2*N-2][N-1:2*N-2]-Y[N-1:2*N-2][N-1:2*N-2])*dz)**2)
'''
	
#**********************************************************************************************************
#Question-4

Q_ij = zeros([2*N-2,2*N-2],dtype=complex)                      #Defining Q_ij & Q_B
Q_B = zeros([2*N-2,1],dtype=complex)

                                                               #Finding Q_ij & Q_B
for i in range(0,2*N-2):
	for j in range(0,2*N-2):
		Q_ij[i][j]=(a/mu0)*((1/(R_ij[i][j]**2))+((1j*k)/R_ij[i][j]))*P_ij[i][j]

Q_B[0:2*N-2]= (a/mu0)*((1/(R_iN[0:2*N-2]**2))+((1j*k)/R_iN[0:2*N-2]))*P_B[0:2*N-2]

#**********************************************************************************************************
#Question-5

J=dot(linalg.inv(M-Q_ij),Im*Q_B)                               #Finding J from M, Q_ij, Q_B

J_final=zeros([2*N+1,1],dtype=complex)                         #Finding J_final from J which has the first, middle and last term too
J_final[0][0]=0
J_final[N][0]=Im
J_final[2*N][0]=0

J_final[1:N]= J[0:N-1]
J_final[N+1:2*N]= J[N-1:]


print('Current obtained through Standard analysis method')
print(I)
print('Current obtained through Magnetic vector potential method')
print(J_final)
title('Comparison of I vs. z')
plot(z,J_final,'g',label = r'Plot from Standard analysis method')
plot(z,I,'r',label = r'Plot from Magnetic vector potential method')
legend()
show()

#Thank you!


