'''
EE2703-Assignment-6
The Laplace Transform
Syam SriBalaji T
EE20B136
17/03/22
'''

from pylab import*
from numpy import*
from math import*
from scipy.signal import*

#***********************************************************************************
#Question-1
#Time response of Spring

Alpha_ = 2.25

Omega1 = 0.5
Alpha1 = 1.5
nume1 = poly1d([1,Omega1])
deno1 = poly1d([1,2*Omega1,Omega1**2+Alpha1**2])
deno1 = polymul(deno1,[1,0,Alpha_])

H1 = lti(nume1,deno1)
t1 = linspace(0,50,1000)
solve1 = impulse(H1,T = t1)
plot(solve1[0],solve1[1])
xlabel('time')
ylabel('F(s)')
title('Spring\'s Time response')
show()

#***********************************************************************************
#Question-2
#Time response of Spring with smaller decay

Omega2 = 0.05
Alpha2 = 1.5
nume2 = poly1d([1,Omega2])
deno2 = poly1d([1,2*Omega2,Omega2**2+Alpha2**2])
deno2 = polymul([1,0,Alpha_],deno2)

H2 = lti(nume2,deno2)
solve2 = impulse(H2,T = t1)
plot(solve2[0],solve2[1],'r')
xlabel('time')
ylabel('F(s)')
title('Spring\'s Time response with smaller decay')
show()

#***********************************************************************************
#Question-3
#LTI response for different frequencies for the applied force

Alpha_ = 2.25
i = 0
for Alpha3 in arange(1.4,1.6,0.05):
	Omega3 = 0.05	
	H3 = lti([1],[1,0,Alpha_])
	t3 = linspace(0,100,1000)
	
	x1 = list(cos(Alpha3*a) for a in t3)
	x2 = list(exp(-1*Omega3*b)*heaviside(b,0.5) for b in t3)
	x_ = multiply(x1,x2)

	f,x,_ = lsim(H3,x_,t3)
	i+=1
	subplot(3,2,i)
	title('frequency = %.2f'%(Alpha3))
	plot(t3,x,'-r')
	

suptitle('LTI response with varying frequency for a force given')
show()

#***********************************************************************************
#Question-4
#Time evolution of Coupled Spring problem

t4 = linspace(0,20,1000)
H_x = lti([1,0,2],[1,0,3,0])
solve_x = impulse(H_x,T = t4)
plot(solve_x[0],solve_x[1],label = 'x')

H_y = lti([2],[1,0,3,0])
solve_y = impulse(H_y,T = t4)
plot(solve_y[0],solve_y[1],label = 'y')

legend()
xlabel('time')
ylabel('x & y')
title('Solution for x(t) & y(t) for Coupled Spring')
show()

#***********************************************************************************
#Question-5
#Steady state Transfer function of Two-port network

H5 = lti([1e6],[1e-6,100,1e6])
w,S,phi = H5.bode()

subplot(2,1,1)
semilogx(w,S,'b')
ylabel(r'$|H(s)|$')

subplot(2,1,2)
semilogx(w,phi,'b')
ylabel(r'$\angle(H(s))$')
suptitle('Magnitude and Phase plot of H(s)')
show()

#***********************************************************************************
#Question-6
#Two-port network with a Input signal

#For 0 < t < 30µs

t6 = linspace(0,30e-6,100)
x6_1 = list(cos(1000*a)-cos(1e6*a) for a in t6)
x6_2 = list(heaviside(b,0.5) for b in t6)
vi_6 = multiply(x6_1,x6_2)
_,y1,svec = lsim(H5,vi_6,t6)
plot(t6,y1,'-g')
xlabel('t')
ylabel('$v_{o}(t)$')
title('Output voltage $v_{o}(t)$ for  0 < t < 30µs')
show()

#For 0 < t < 10ms

t7 = linspace(0,10*1e-3,100000)
x7_1 = list(cos(1000*a)-cos(1e6*a) for a in t7)
x7_2 = list(heaviside(b,0.5) for b in t7)
vi_7 = multiply(x7_1,x7_2)
_,y2,svec = lsim(H5,vi_7,t7)
plot(t7,y2,'-g')
xlabel('t')
ylabel('$v_{o}(t)$')
title('Output voltage $v_{o}(t)$ for  0 < t < 10ms')
show()











