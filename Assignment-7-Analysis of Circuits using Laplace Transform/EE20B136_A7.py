'''
EE2703-Assignment-7
Circuit analysis using SymPy
Syam SriBalaji T
EE20B136
31/03/22
'''

from sympy import *
from pylab import *
from scipy.signal import *
from numpy import *

def SymPy_to_LTI(sympy_expr,s = symbols('s')):
	nume, deno = simplify(sympy_expr).as_numer_denom()
	c_nume = Poly(nume,s).all_coeffs()
	c_deno = Poly(deno,s).all_coeffs()
	c_nume = [float(i) for i in c_nume]
	c_deno = [float(i) for i in c_deno]
	return lti(c_nume, c_deno)

#****************************************************************************************
#Low pass filter

s = symbols('s')
def Lowpass(R1,R2,C1,C2,G,Vi):
	A1 = Matrix([[0,0,1,-1/G],[-1/(1+s*R2*C2),1,0,0],[0,-G,G,1],[-1/R1-1/R2-s*C1,1/R2,0,s*C1]])
	b1 = Matrix([0,0,0,-Vi/R1])
	V1 = A1.inv()*b1
	return (A1,b1,V1)

#****************************************************************************************
#Magnitude response of the Low pass filter
#Question:1

A1,b1,V1 = Lowpass(10000,10000,1e-9,1e-9,1.586,1)
Vo_LP_1 = V1[3]
Omega1 = logspace(0,8,1000)
S1 = 1j*Omega1
H_func1 = lambdify(s,Vo_LP_1,'numpy')
Values_1 = H_func1(S1)

loglog(Omega1,abs(Values_1),'r')
xlabel('ω')
ylabel(r'$|H(jω)|$')
title('Magnitude response of LPF ( $|H(jω)|$ vs ω ) in loglog')
grid(True)
show()

#****************************************************************************************
#Mixed frequency sinusoid response of Low pass filter
#Question:2

H1 = SymPy_to_LTI(Vo_LP_1)
t1 = linspace(0,0.01,100000)
Vi = multiply((sin(2000*pi*t1)+cos(2000000*pi*t1)),heaviside(t1,0.5))
Vo_LP_3 = lsim(H1,Vi,T = t1)

plot(Vo_LP_3[0],Vi,label = r'$V_{in}$')
plot(Vo_LP_3[0],Vo_LP_3[1],'r',label = r'$V_{o}$')
xlabel('t')
ylabel('V')
title('Output response of Mixed frequency sinusoidal inputs')
grid(True)
legend()
show()

#****************************************************************************************
#Step response of the Low pass filter
#Question:1

t2 = linspace(0,0.001,1000)
Vo_LP_2 = step(H1,T = t2)

plot(Vo_LP_2[0],Vo_LP_2[1],'g')
xlabel('t')
ylabel('Vo')
title('Step response of the Low pass filter')
grid(True)
show()

#****************************************************************************************
#High pass filter

def Highpass(R1,R3,C1,C2,G,Vi):
	A2 = Matrix([[0,0,1,-1/G],[-1/(1+1/(s*R3*C2)),1,0,0],[0,-G,G,1],[-s*C1-s*C2-1/R1,s*C2,0,1/R1]])
	b2 = Matrix([0,0,0,-Vi*s*C1])
	V2 = A2.inv()*b2
	return (A2,b2,V2)

#****************************************************************************************
#Magnitude response of the High pass filter
#Question:3

A2,b2,V2 = Highpass(10000,10000,1e-9,1e-9,1.586,1)
Vo_HP_1 = V2[3]
Omega2 = logspace(0,8,1000)
S2 = 1j*Omega2
H_func2 = lambdify(s,Vo_HP_1,'numpy')
V2 = H_func2(S2)

loglog(Omega2,abs(V2),'r')
xlabel('ω')
ylabel(r'$|H(jω)|$')
title('Magnitude response of HPF ( $|H(jω)|$ vs ω ) in loglog')
grid(True)
show()

#****************************************************************************************
#Response of Damped sinusoidal with input of 1 Hz in HPF
#Question:4-(i)

H2 = SymPy_to_LTI(Vo_HP_1)
t3 = linspace(0,10,1000)
Vi = multiply(multiply(exp(-0.5*t3),sin(2*pi*t3)),heaviside(t3,0.5))
Vo_HP_2 = lsim(H2,Vi,T = t3)

plot(Vo_HP_2[0],Vi,'b',label = r'$V_{in}$')
plot(Vo_HP_2[0],Vo_HP_2[1],'r',label = r'$V_{o}$')
xlabel('t')
ylabel('V')
title('Response of Damped sinusoidal with input of 1 Hz in HPF')
grid(True)
legend()
show()

#****************************************************************************************
#Response of Damped sinusoidal with input of 2x10^5 Hz in HPF
#Question:4-(ii)

t4 = linspace(0,0.0001,10000)
Vi = multiply(multiply(exp(-0.5*t4),sin(2*pi*200000*t4)),heaviside(t4,0.5))
Vo_HP_3 = lsim(H2,Vi,T = t4)

plot(Vo_HP_3[0],Vi,'b',label = r'$V_{in}$')
plot(Vo_HP_3[0],Vo_HP_3[1],'r',label = r'$V_{o}$')
xlabel('t')
ylabel('V')
title('Response of Damped sinusoidal with input of $2x10^5$ Hz in HPF')
grid(True)
legend()
show()

#****************************************************************************************
#Step response for High pass filter
#Question:5

t5 = linspace(0,0.001,1000)
Vo_HP_4 = step(H2,T = t5)

plot(Vo_HP_4[0],Vo_HP_4[1],'g')
xlabel('t')
ylabel('$V_{o}$')
title('Step response of the High pass filter')
grid(True)
show()

#Thank you!


