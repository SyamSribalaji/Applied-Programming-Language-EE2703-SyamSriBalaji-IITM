'''
EE2703-Assignment-8
Digital Fourier Transform
Syam SriBalaji T
EE20B136
14/04/22
'''

from sympy import *
from pylab import *
from numpy import *

#****************************************************************************************
#Spectrum of sin^3(t) and cos^3(t)

t=linspace(-4*pi, 4*pi, 513);t=t[:-1]
y1=sin(t)**3
y2=cos(t)**3
Y1=fftshift(fft.fft(y1))/512.0
Y2=fftshift(fft.fft(y2))/512.0
w=linspace(-64, 64, 513);w=w[:-1]

subplot(2,1,1)
plot(w, abs(Y1), lw=2)
xlim([-15,15])
ylabel(r"$|X|$", size=16)
title(r"Spectrum of $sin^3(t)$")
grid(True)
subplot(2,1,2)
plot(w, angle(Y1), 'ko', lw=2, mfc='r', markeredgewidth=1.0)
xlim([-15,15])
ii1 = where(abs(Y1)>1e-3)
plot(w[ii1], angle(Y1[ii1]), 'ko', lw=2, mfc='g', markeredgewidth=1.0)
ylabel(r"$∠X$", size=16)
xlabel(r"$\omega$", size=16)
grid(True)
show()

subplot(2,1,1)
plot(w, abs(Y2), lw=2)
xlim([-15,15])
ylabel(r"$|X|$", size=16)
title(r"Spectrum of $cos^3(t)$")
grid(True)
subplot(2,1,2)
plot(w, angle(Y2), 'ko', lw=2, mfc='r', markeredgewidth=1.0)
xlim([-15,15])
ii2 = where(abs(Y2)>1e-3)
plot(w[ii2], angle(Y2[ii2]), 'ko', lw=2, mfc='g', markeredgewidth=1.0)
ylabel(r"∠X", size=16)
xlabel(r"$\omega$", size=16)
grid(True)
show()

#****************************************************************************************
#Spectrum of cos(20t +5 cos(t))

y3=cos(20*t+5*cos(t))
Y3=fftshift(fft.fft(y3))/512.0

subplot(2,1,1)
plot(w, abs(Y3), lw=2)
xlim([-30,30])
ylabel(r"$|X|$", size=16)
title(r"Spectrum of $cos(20t+5cos(t))$")
grid(True)
subplot(2,1,2)
xlim([-30,30])
ii3 = where(abs(Y3)>1e-3)
plot(w[ii3], angle(Y3[ii3]), 'ko', lw=2, mfc='g', markeredgewidth=1.0)
ylabel(r"$∠X$", size=16)
xlabel(r"$\omega$", size=16)
grid(True)
show()

#****************************************************************************************
#Spectrum and Maximum error of exp(-t^2/2)

def Gauss(r,N):
	t2=linspace(-r*pi, r*pi, N+1);t2=t2[:-1]
	y4=exp(-1*(t2**2)/2)
	Y4= fftshift(fft.fft(y4))/N
	w=linspace(-32, 32, N+1);w=w[:-1]
	
	Y4= fftshift(abs(fft.fft(y4)))/N
	Y4= Y4*sqrt(2*pi)/max(Y4)
	Y4_=exp(-w**2/2)*sqrt(2*pi)

	subplot(2,1,1)
	annotate("Maximum error = %e" %(abs(Y4-Y4_).max()), xy =(3.3, 1),xytext =(2.5, 2.25))
	plot(w, abs(Y4), lw=2)
	xlim([-10,10])
	ylabel(r"$|X|$", size=16)
	title(r"Spectrum of $e^{(-t^2/2)}$ in (-%d$\pi$,%d$\pi$) & N=%d"%(r,r,N))
	grid(True)

	subplot(2,1,2)
	plot(w, angle(Y4), 'ko', lw=2, mfc='r', markeredgewidth=1.0)
	ii4 = where(abs(Y4)>1e-3)
	plot(w[ii4], angle(Y4[ii4]), 'ko', lw=2, mfc='g', markeredgewidth=1.0)
	xlim([-10,10])
	ylabel(r"$∠X$", size=16)
	xlabel(r"$\omega$", size=16)
	grid(True)
	show()
	print("Maximum error is {}".format(abs(Y4-Y4_).max()))

'''
Gauss(2,256)
Gauss(2,512)
Gauss(2,1024)
Gauss(4,256)
Gauss(4,512)
Gauss(4,1024)
Gauss(8,256)
Gauss(8,512)
Gauss(8,1024)
'''
print("DFT of exp(-t^2/2) has maximum accuracy for range (-4π,4π) and N=256")
Gauss(4,256)

#Thank you!

