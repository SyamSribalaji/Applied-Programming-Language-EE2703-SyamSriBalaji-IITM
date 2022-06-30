'''
EE2703-Assignment-4
Fourier Approximations
Syam SriBalaji T
EE20B136
25/02/22
'''

from pylab import*
from scipy.special import*
from numpy import*
from math import*
from sys import*
from scipy import integrate

PI=pi
def E(e):
	return exp(e)
def CC(c):
	return cos(cos(c))

#***********************************************************************************
#Question-1

#Plotting Figure-1
#e^x Function from (-2π,4π)
N=1000
x_coeff1=linspace(-2*PI,4*PI,N)
y_coeff1=list(0 for c in range(0,N))
for j in range(0,N):
	y_coeff1[j]=E(x_coeff1[j])
plot(x_coeff1,y_coeff1)
xlabel('x-axis')
ylabel('y-axis')
grid(True)
title('Function 1: e^x Graph')
show()

#Plotting Figure-2
#cos(cos(x)) Function from (-2π,4π)
y_coeff2=list(0 for c in range(0,N))
for j in range(0,N):
	y_coeff2[j]=CC(x_coeff1[j])
plot(x_coeff1,y_coeff2)
xlabel('x-axis')
ylabel('y-axis')
grid(True)
title('Function 2: Cos(Cos(x))')
show()

#***********************************************************************************
#Question-2

#Finding Coefficients of e^x function in Integration method

f_func1 = lambda x: E(x)
g_func1 = integrate.quad(f_func1, 0, 2*PI)
Ao_1=(1/(2*PI))*g_func1[0]
N1=25

A_coeff1=zeros([N1,1])
for i in range(0,N1):
	k=i+1
	f_func2= lambda x,k: E(x)*cos(k*x)
	g_func2 = integrate.quad(f_func2, 0, 2*PI,args=(k))
	A_coeff1[i][0]=(1/PI)*g_func2[0]

B_coeff1=zeros([N1,1])
for i in range(0,N1):
	k=i+1
	f_func3= lambda x,k: E(x)*sin(k*x)
	g_func3 = integrate.quad(f_func3, 0, 2*PI,args=(k))
	B_coeff1[i][0]=(1/PI)*g_func3[0]

#Finding Coefficients of cos(cos(x)) function in Integration method

f_func4 = lambda x: CC(x)
g_func4 = integrate.quad(f_func4, 0, 2*PI)
Ao_2=(1/(2*PI))*g_func4[0]

A_coeff2=zeros([N1,1])
for i in range(0,N1):
	k=i+1
	f_func5= lambda x,k: CC(x)*cos(k*x)
	g_func5 = integrate.quad(f_func5, 0, 2*PI,args=(k))
	A_coeff2[i][0]=(1/PI)*g_func5[0]
	
B_coeff2=zeros([N1,1])
for i in range(0,N1):
	k=i+1
	f_func6= lambda x,k: CC(x)*sin(k*x)
	g_func6 = integrate.quad(f_func6, 0, 2*PI,args=(k))
	B_coeff2[i][0]=(1/PI)*g_func6[0]

#***********************************************************************************
#Question-3

#Creating the Coefficient vector for e^x function by integration method
#coeff_vect1 is created here

coeff_vect1 = zeros(2*N1+1)
coeff_vect1[0]=Ao_1
t=0
for s in range(1,2*N1,2):
	coeff_vect1[s]=A_coeff1[t][0]
	coeff_vect1[s+1]=B_coeff1[t][0]
	t+=1

#Creating the Coefficient vector for cos(cos(x)) function by integration method
#coeff_vect2 is created here

coeff_vect2 = zeros(2*N1+1)
coeff_vect2[0]=Ao_2
t=0
for s in range(1,2*N1,2):
	coeff_vect2[s]=A_coeff2[t][0]
	coeff_vect2[s+1]=B_coeff2[t][0]
	t+=1

#Plotting Figure-3

plot(0,abs(Ao_1),'ro')
semilogy(range(1,N1+1),abs(coeff_vect1[1::2]),'ro',label='An')
semilogy(range(1,N1+1),abs(coeff_vect1[2::2]),'bo',label='Bn')
xlabel('n values')
ylabel('A & B for e^x function')
grid(True)
title('semilog A & B for e^x function')
legend()
show()

#Plotting Figure-4

plot(0,abs(Ao_1),'ro')
loglog(range(1,N1+1),abs(coeff_vect1[1::2]),'ro',label='An')
loglog(range(1,N1+1),abs(coeff_vect1[2::2]),'bo',label='Bn')
xlabel('n values')
ylabel('A & B for e^x function')
grid(True)
title('loglog A & B for e^x function')
legend()
show()

#Plotting Figure-5

plot(0,abs(Ao_2),'ro')
semilogy(range(1,N1+1),abs(coeff_vect2[1::2]),'ro',label='An')
semilogy(range(1,N1+1),abs(coeff_vect2[2::2]),'bo',label='Bn')
xlabel('n values')
ylabel('A & B for cos(cos(x)) function')
grid(True)
title('semilog A & B for cos(cos(x)) function')
legend()
show()

#Plotting Figure-6

plot(0,abs(Ao_2),'ro')
loglog(range(1,N1+1),abs(coeff_vect2[1::2]),'ro',label='An')
loglog(range(1,N1+1),abs(coeff_vect2[2::2]),'bo',label='Bn')
xlabel('n values')
ylabel('A & B for cos(cos(x)) function')
grid(True)
title('loglog A & B for cos(cos(x)) function')
legend()
show()

#***********************************************************************************
#Question-4

#Finding Coefficient Vector for e^x funciton in Least squares approch method
N2=300
x=linspace(0,2*pi,N2+1)
x=x[:-1]

b1=zeros([N2,1])
for n in range(0,N2):
	b1[n,0]=E(x[n])

A1=zeros((N2,2*N1+1))
A1[:,0]=1
for k in range(1,N1+1):
	for m in range(0,N2):
		A1[m,(2*k)-1]=cos(k*x[m])
		A1[m,2*k]=sin(k*x[m])
c1=lstsq(A1,b1,rcond=None)[0]

#***********************************************************************************
#Finding Coefficient Vector for cos(cos(x)) function in Least squares approch method

b2=zeros([N2,1])
for n in range(0,N2):
	b2[n,0]=CC(x[n])

A2=zeros((N2,2*N1+1))
A2[:,0]=1
for k in range(1,N1+1):	
	for m in range(0,N2):
		A2[m,(2*k)-1]=cos(k*x[m])
		A2[m,2*k]=sin(k*x[m])
c2=lstsq(A2,b2,rcond=None)[0]

#***********************************************************************************
#Plotting Figure-7
#Comparison of Semilog A & B coefficients in e^x function in both methods

plot(0,abs(Ao_1),'ro')
semilogy(range(1,N1+1),abs(coeff_vect1[1::2]),'ro',label='An in Integration method')
semilogy(range(1,N1+1),abs(coeff_vect1[2::2]),'bo',label='Bn in Integration method')
plot(0,abs(c1[0]),'go')
semilogy(range(1,N1+1),abs(c1[1::2]),'go',label='An in Least squares approch method')
semilogy(range(1,N1+1),abs(c1[2::2]),'yo',label='Bn in Least squares approch method')
xlabel('n values')
ylabel('A & B for e^x function')
grid(True)
title('Comparison of semilog A & B in e^x function in both methods')
legend()
show()

#***********************************************************************************
#Plotting Figure-8
#Comparison of loglog A & B coefficients in e^x function in both methods

plot(0,abs(Ao_1),'ro')
loglog(range(1,N1+1),abs(coeff_vect1[1::2]),'ro',label='A in Integration method')
loglog(range(1,N1+1),abs(coeff_vect1[2::2]),'bo',label='B in Integration method')
plot(0,abs(c1[0]),'go')
loglog(range(1,N1+1),abs(c1[1::2]),'go',label='An in Least squares approch method')
loglog(range(1,N1+1),abs(c1[2::2]),'yo',label='Bn in Least squares approch method')
xlabel('n values')
ylabel('A & B for e^x function')
title('Comparison of loglog A & B in e^x function in both methods')
legend()
show()

#***********************************************************************************
#Plotting Figure-9
#Comparison of Semilog A & B coefficients in cos(cos(x)) function in both methods

plot(0,abs(Ao_2),'ro')
semilogy(range(1,N1+1),abs(coeff_vect2[1::2]),'ro',label='An in Integration method')
semilogy(range(1,N1+1),abs(coeff_vect2[2::2]),'bo',label='Bn in Integration method')
plot(0,abs(c2[0]),'go')
semilogy(range(1,N1+1),abs(c2[1::2]),'go',label='An in Least squares approch method')
semilogy(range(1,N1+1),abs(c2[2::2]),'yo',label='Bn in Least squares approch method')
xlabel('n values')
ylabel('A & B for cos(cos(x)) function')
grid(True)
title('Comparison of semilog A & B in cos(cos(x)) function in both methods')
legend()
show()

#***********************************************************************************
#Plotting Figure-10
#Comparison of loglog A & B coefficients in cos(cos(x)) function in both methods

plot(0,abs(Ao_2),'ro')
loglog(range(1,N1+1),abs(coeff_vect2[1::2]),'ro',label='An in Integration method')
loglog(range(1,N1+1),abs(coeff_vect2[2::2]),'bo',label='Bn in Integration method')
plot(0,abs(c2[0]),'go')
loglog(range(1,N1+1),abs(c2[1::2]),'go',label='An in Least squares approch method')
loglog(range(1,N1+1),abs(c2[2::2]),'yo',label='Bn in Least squares approch method')
xlabel('n values')
ylabel('A & B for cos(cos(x)) function')
title('Comparison of loglog A & B in cos(cos(x)) function in both methods')
legend()
show()

#***********************************************************************************
#Question-6
#Finding Absolute difference between the two sets of coefficients and its largest deviation

CR1=c1.reshape((2*N1+1,))
CR2=c2.reshape((2*N1+1,))
Dev_E=coeff_vect1- CR1
Dev_CC=coeff_vect2- CR2
Max_Dev_E= max(Dev_E)
Max_Dev_CC= max(Dev_CC)
print('\n\nMaximum deviation of coefficients in e^x function is %e'%(Max_Dev_E))
print('Maximum deviation of coefficients in cos(cos(x)) function is %e'%(Max_Dev_CC))

#***********************************************************************************
#Question-7
#Plotting the e^x function values got from Least squares approch method

Approx_LSA_E= dot(A1,c1)
Approx_LSA_CC= dot(A2,c2)

#Plotting Figure-11

plot(x,Approx_LSA_E,'go',label='e^x graph obtained by Matrix multiplication')
plot(x,b1,'r',label='True e^x function')
xlabel('n values')
ylabel('e^x function values')
grid(True)
title('Comparison of Approximate e^x graph obtained with True function')
legend()
show()

#Plotting Figure-12

plot(x,Approx_LSA_CC,'go',label='cos(cos(x)) graph obtained by Matrix multiplication')
plot(x,b2,'r',label='True cos(cos(x)) function')
xlabel('n values')
ylabel('cos(cos(x)) function values')
grid(True)
title('Comparison of Approximate cos(cos(x)) graph obtained with True function')
legend()
show()

#***********************************************************************************
#Question-1
#Making Periodic extension & e^x by Least Square Approch and Integration method
#Plotting Figure-13

Approx_LSA_E= dot(A1,c1)
Approx_I_E= dot(A1,coeff_vect1)

x_coeff_PE= linspace(-2*PI,4*PI,3*N2)
x1_coeff_PE = linspace(0,2*PI,N2)
x2_coeff_PE = linspace(2*PI,4*PI,N2)
x3_coeff_PE = linspace(-2*PI,0,N2)
E_y_PE = list(E(c) for c in x1_coeff_PE)
E1_y_PE=list(E(d) for d in x_coeff_PE)

plot(x1_coeff_PE,Approx_I_E,'bo',label='e^x graph by Integration method')
plot(x2_coeff_PE,Approx_I_E,'bo')
plot(x3_coeff_PE,Approx_I_E,'bo')
plot(x1_coeff_PE,Approx_LSA_E,'go',label='e^x graph by Least Square Approch method')
plot(x2_coeff_PE,Approx_LSA_E,'go')
plot(x3_coeff_PE,Approx_LSA_E,'go')

semilogy(x_coeff_PE,E1_y_PE,'b',label='Correct value')
semilogy(x1_coeff_PE,E_y_PE,'r',label='Periodic Extension')
semilogy(x2_coeff_PE,E_y_PE,'r')
semilogy(x3_coeff_PE,E_y_PE,'r')

grid(True)
ylabel('e^x')
xlabel('x')
legend()
show()

#Thank you!














	











