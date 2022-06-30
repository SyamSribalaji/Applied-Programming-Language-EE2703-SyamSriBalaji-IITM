'''
EE2703-Assignment-9
Spectra of Non-Periodic signals
Syam SriBalaji T
EE20B136
29/04/22
'''

from sympy import *
from numpy import *
from pylab import *

#****************************************************************************************
#Spectrum of sin(√2t)

t = linspace(-pi, pi, 65);t = t[:-1]
dt = t[1]-t[0];fmax = 1/dt
y = sin(sqrt(2)*t)
y[0] = 0
y = fftshift(y)
Y = fftshift(fft(y))/64.0
w = linspace(-pi*fmax, pi*fmax, 65);w = w[:-1]
subplot(2, 1, 1)
plot(w, abs(Y), lw = 2)
xlim([-10, 10])
ylabel(r"$|Y|$", size = 16)
title(r"Spectrum of $\sin\left(\sqrt{2}t\right)$")
grid(True)
subplot(2, 1, 2)
plot(w, angle(Y), 'ko', lw = 2, mfc = 'r', markeredgewidth = 1.0)
xlim([-10, 10])
ylabel(r"Phase of $Y$", size = 16)
xlabel(r"$\omega$", size = 16)
grid(True)
show()

#****************************************************************************************
#Plot of sin(√2t)

t1 = linspace(-pi, pi, 65);t1 = t1[:-1]
t2 = linspace(-3*pi, -pi, 65);t2 = t2[:-1]
t3 = linspace(pi, 3*pi, 65);t3 = t3[:-1]

plot(t1, sin(sqrt(2)*t1), 'b', lw = 2)
plot(t2, sin(sqrt(2)*t2), 'r', lw = 2)
plot(t3, sin(sqrt(2)*t3), 'r', lw = 2)
ylabel(r"$y$", size = 16)
xlabel(r"$t$", size = 16)
title(r"Plot of $\sin\left(\sqrt{2}t\right)$")
grid(True)
show()

#****************************************************************************************
#Discrete plot of sin(√2t)

t1 = linspace(-pi, pi, 65);t1 = t1[:-1]
t2 = linspace(-3*pi, -pi, 65);t2 = t2[:-1]
t3 = linspace(pi, 3*pi, 65);t3 = t3[:-1]
y = sin(sqrt(2)*t1)
plot(t1, y, 'ko', lw = 2, mfc = 'b', markeredgewidth = 1.0)
plot(t2, y, 'ko', lw = 2, mfc = 'r', markeredgewidth = 1.0)
plot(t3, y, 'ko', lw = 2, mfc = 'r', markeredgewidth = 1.0)
ylabel(r"$y$", size = 16)
xlabel(r"$t$", size = 16)
title(r"$\sin\left(\sqrt{2}t\right)$ with $t$ wrapping every $2\pi$ ")
grid(True)
show()

#****************************************************************************************
#Spectrum of Periodic digital ramp

t = linspace(-pi, pi, 65);t = t[:-1]
dt = t[1]-t[0];fmax = 1/dt
y = t
y[0] = 0
y = fftshift(y)
Y = fftshift(fft(y))/64.0
w = linspace(-pi*fmax, pi*fmax, 65);w = w[:-1]
semilogx(abs(w), 20*log10(abs(Y)), lw = 2)
xlim([1, 10])
ylim([-20, 0])
ylabel(r"$|Y|$ (dB)", size = 16)
title(r"Spectrum of a Periodic digital ramp")
xlabel(r"$\omega$", size = 16)
grid(True)
show()

#****************************************************************************************
#sin(√2t)×ω(t) with time gaps of 2π

t1 = linspace(-pi, pi, 65);t1 = t1[:-1]
t2 = linspace(-3*pi, -pi, 65);t2 = t2[:-1]
t3 = linspace(pi, 3*pi, 65);t3 = t3[:-1]
n = arange(64)
wnd = fftshift(0.54+0.46*cos(2*pi*n/63))
y = sin(sqrt(2)*t1)*wnd

plot(t1, y, 'ko', lw = 2, mfc = 'b', markeredgewidth = 1.0)
plot(t2, y, 'ko', lw = 2, mfc = 'r', markeredgewidth = 1.0)
plot(t3, y, 'ko', lw = 2, mfc = 'r', markeredgewidth = 1.0)
ylabel(r"$y$", size = 16)
xlabel(r"$t$", size = 16)
title(r"$\sin\left(\sqrt{2}t\right)\times ω(t)$ with $t$ wrapping every $2\pi$ ")
grid(True)
show()

#****************************************************************************************
#Spectrum of sin(√2t)×ω(t) with time period 2π

t = linspace(-pi, pi, 65);t = t[:-1]
dt = t[1]-t[0];fmax = 1/dt
n = arange(64)
wnd = fftshift(0.54+0.46*cos(2*pi*n/63))
y = sin(sqrt(2)*t)*wnd
y[0] = 0
y = fftshift(y)
Y = fftshift(fft(y))/64.0
w = linspace(-pi*fmax, pi*fmax, 65);w = w[:-1]

subplot(2, 1, 1)
plot(w, abs(Y), lw = 2)
xlim([-8, 8])
ylabel(r"$|Y|$", size = 16)
title(r"Spectrum of $\sin\left(\sqrt{2}t\right)\times ω(t)$ $with$ time period $2\pi$")
grid(True)
subplot(2, 1, 2)
plot(w, angle(Y), 'ko', lw = 2, mfc = 'r', markeredgewidth = 1.0)
xlim([-8, 8])
ylabel(r"Phase of $Y$", size = 16)
xlabel(r"$\omega$", size = 16)
grid(True)
show()

#****************************************************************************************
#Spectrum of sin(√2t)×ω(t) with time period 8π

t = linspace(-4*pi, 4*pi, 257);t = t[:-1]
dt = t[1]-t[0];fmax = 1/dt
n = arange(256)
wnd = fftshift(0.54+0.46*cos(2*pi*n/256))
y = sin(sqrt(2)*t)
y = y*wnd
y[0] = 0
y = fftshift(y)
Y = fftshift(fft(y))/256.0
w = linspace(-pi*fmax, pi*fmax, 257);w = w[:-1]

subplot(2, 1, 1)
plot(w, abs(Y))
xlim([-4, 4])
ylabel(r"$|Y|$", size = 16)
title(r"Improved Spectrum of $\sin\left(\sqrt{2}t\right)\times ω(t)$ $with$ time period $8\pi$")
grid(True)
subplot(2, 1, 2)
plot(w, angle(Y), 'ko', lw = 2, mfc = 'r', markeredgewidth = 1.0)
xlim([-4, 4])
ylabel(r"Phase of $Y$", size = 16)
xlabel(r"$\omega$", size = 16)
grid(True)
show()


#****************************************************************************************
#Question-2
#FFT of cos^3(0.86t) without Windowing

t = linspace(-4*pi, 4*pi, 257)[:-1]
dt = t[1]-t[0];
fmax = 1/dt
y = (cos(0.86*t))**3
y[0] = 0
y = fftshift(y)
Y = fftshift(fft(y))/256
w = linspace(-fmax*pi, fmax*pi, 257)[:-1]
    
Mag = abs(Y)
Ph = angle(Y)
    
subplot(2, 1, 1)
plot(w, Mag, lw = 2)
xlim([-3, 3])
ylabel(r"$|Y|$", size = 16)
title(r"Spectrum of $cos^3(ω_0t)$ without Windowing")
grid(True)
subplot(2, 1, 2)
Ph[where(Mag<3e-3)] = 0
plot(w, Ph, 'ko', lw = 2, mfc = 'r', markeredgewidth = 1.0)
xlim([-3, 3])
ylabel(r"Phase of $Y$", size = 16)
xlabel(r"$\omega$", size = 16)
grid(True)
show()

#****************************************************************************************
#FFT of cos^3(0.86t) with Windowing

t = linspace(-4*pi, 4*pi, 257)[:-1]
dt = t[1]-t[0];
fmax = 1/dt
y = (cos(0.86*t))**3

m = arange(256)
wnd = fftshift(0.54+0.46*cos(2*pi*m/256))
y1 = y*wnd
    
y1[0] = 0
y1 = fftshift(y1)
Y1 = fftshift(fft(y1))/256
w = linspace(-fmax*pi, fmax*pi, 257)[:-1]
    
Mag = abs(Y1)
Ph = angle(Y1)
    
subplot(2, 1, 1)
plot(w, Mag, lw = 2)
xlim([-3, 3])
ylabel(r"$|Y|$", size = 16)
title(r"Spectrum of $cos^3(ω_0t)$ with Windowing")
grid(True)
subplot(2, 1, 2)
Ph[where(Mag<3e-3)] = 0
plot(w, Ph, 'ko', lw = 2, mfc = 'r', markeredgewidth = 1.0)
xlim([-3, 3])
ylabel(r"Phase of $Y$", size = 16)
xlabel(r"$\omega$", size = 16)
grid(True)
show()


#****************************************************************************************
#Question-3
#FFT of cos(1.5t+0.5)

t = linspace(-pi, pi, 129)[:-1]
dt = t[1]-t[0];
fmax = 1/dt
y = cos(1.5*t +0.5)

m = arange(128)
wnd = fftshift(0.54+0.46*cos(2*pi*m/128))
y = y*wnd
    
y[0] = 0
y = fftshift(y)
Y = fftshift(fft(y))/float(128)
w = linspace(-pi*fmax, pi*fmax, 129)[:-1]
    
Mag = abs(Y)
Ph = angle(Y)
subplot(2, 1, 1)
plot(w, Mag, lw = 2)
xlim([-3, 3])
ylabel(r"$|Y|$", size = 16)
title(r"Spectrum of $cos(ω_0t + \delta)$")
grid(True)
subplot(2, 1, 2)
Ph[where(Mag<3e-3)] = 0
plot(w, Ph, 'ko', lw = 2, mfc = 'r', markeredgewidth = 1.0)
xlim([-3, 3])
ylabel(r"Phase of $Y$", size = 16)
xlabel(r"$\omega$", size = 16)
grid(True)
show()

ii = where(w>0)
omega = (sum(abs(Y[ii])**2*w[ii])/sum(abs(Y[ii])**2))
print ("Omega = ", omega)

ii_1 = where(logical_and(abs(Y)>1e-4, w>0))[0]
sort(ii_1)
points = ii_1[1:2]
print ("Delta = ", sum(angle(Y[points]))/len(points))

#****************************************************************************************
#Question-4
#FFT of cos(1.5t+0.5) with noise

t = linspace(-pi, pi, 129)[:-1]
dt = t[1]-t[0];
fmax = 1/dt
y = cos(1.5*t + 0.5) + 0.1*np.random.randn(128)

m = arange(128)
wnd = fftshift(0.54+0.46*cos(2*pi*m/128))
y = y*wnd
    
y[0] = 0
y = fftshift(y)
Y = fftshift(fft(y))/float(128)
w = linspace(-pi*fmax, pi*fmax, 129)[:-1]
    
Mag = abs(Y)
Ph = angle(Y)
subplot(2, 1, 1)
plot(w, Mag, lw = 2)
xlim([-3, 3])
ylabel(r"$|Y|$", size = 16)
title(r"Spectrum of $cos(ω_0t + \delta)$")
grid(True)
subplot(2, 1, 2)
Ph[where(Mag<3e-3)] = 0
plot(w, Ph, 'ko', lw = 2, mfc = 'r', markeredgewidth = 1.0)
xlim([-3, 3])
ylabel(r"Phase of $Y$", size = 16)
xlabel(r"$\omega$", size = 16)
grid(True)
show()

ii = where(w>0)
omega = (sum(abs(Y[ii])**2*w[ii])/sum(abs(Y[ii])**2))
print ("Omega = ", omega)

ii_1 = where(logical_and(abs(Y)>1e-4, w>0))[0]
sort(ii_1)
points = ii_1[1:2]
print ("Delta = ", sum(angle(Y[points]))/len(points))

#****************************************************************************************
#Question-5
#cos(16t(1.5 + t/2π)) with Windowing

t = linspace(-pi, pi, 1025)[:-1]
dt = t[1]-t[0];
fmax = 1/dt
y = cos(16*(1.5 + t/(2*pi))*t)

m = arange(1024)
wnd = fftshift(0.54+0.46*cos(2*pi*m/1024))
y = y*wnd
    
y[0] = 0
y = fftshift(y)
Y = fftshift(fft(y))/float(1024)
w = linspace(-pi*fmax, pi*fmax, 1025)[:-1]
    
Mag = abs(Y)
Ph = angle(Y)
subplot(2, 1, 1)
plot(w, Mag, lw = 2)
xlim([-60, 60])
ylabel(r"$|Y|$", size = 16)
title(r"Spectrum of Chirp function with Windowing")
grid(True)
subplot(2, 1, 2)
Ph[where(Mag<3e-3)] = 0
plot(w, Ph, 'ko', lw = 2, mfc = 'r', markeredgewidth = 1.0)
xlim([-60, 60])
ylabel(r"Phase of $Y$", size = 16)
xlabel(r"$\omega$", size = 16)
grid(True)
show()

ii = where(w>0)
omega = (sum(abs(Y[ii])**2*w[ii])/sum(abs(Y[ii])**2))
print ("Omega = ", omega)

ii_1 = where(logical_and(abs(Y)>1e-4, w>0))[0]
sort(ii_1)
points = ii_1[1:2]
print ("Delta = ", sum(angle(Y[points]))/len(points))

#****************************************************************************************
#cos(16t(1.5 + t/2π)) without Windowing

t = linspace(-pi, pi, 1025)[:-1]
dt = t[1]-t[0];
fmax = 1/dt
y = cos(16*(1.5 + t/(2*pi))*t)

y[0] = 0
y = fftshift(y)
Y = fftshift(fft(y))/float(1024)
w = linspace(-pi*fmax, pi*fmax, 1025)[:-1]
    
Mag = abs(Y)
Ph = angle(Y)
subplot(2, 1, 1)
plot(w, Mag, lw = 2)
xlim([-60, 60])
ylabel(r"$|Y|$", size = 16)
title(r"Spectrum of Chirp function without Windowing")
grid(True)
subplot(2, 1, 2)
Ph[where(Mag<3e-3)] = 0
plot(w, Ph, 'ko', lw = 2, mfc = 'r', markeredgewidth = 1.0)
xlim([-60, 60])
ylabel(r"Phase of $Y$", size = 16)
xlabel(r"$\omega$", size = 16)
grid(True)
show()

ii = where(w>0)
omega = (sum(abs(Y[ii])**2*w[ii])/sum(abs(Y[ii])**2))
print ("Omega = ", omega)

ii_1 = where(logical_and(abs(Y)>1e-4, w>0))[0]
sort(ii_1)
points = ii_1[1:2]
print ("Delta = ", sum(angle(Y[points]))/len(points))


#****************************************************************************************
#Question-6
#Chirp function 

t = linspace(-pi, pi, 1025);t = t[:-1]
t_arrays = split(t, 16)

Y_Mag = zeros((16, 64))
Y_Angle = zeros((16, 64))

Omega_list = list(0 for k in range(len(t_arrays)))
Delta_list = list(0 for k in range(len(t_arrays)))

for i in range(len(t_arrays)):
    t = t_arrays[i]
    dt = t[1]-t[0];
    fmax = 1/dt
    y = cos(16*(1.5 + t/(2*pi))*t)
    y[0] = 0
    y = fftshift(y)
    Y = fftshift(fft(y))/float(64)
    w = linspace(-pi*fmax, pi*fmax, 65)[:-1]
    Mag = abs(Y)
    Ph = angle(Y)
    ii = where(w>0)
    omega = (sum(abs(Y[ii])**2*w[ii])/sum(abs(Y[ii])**2))
    Omega_list[i] = omega
    ii_1 = where(logical_and(abs(Y)>1e-4, w>0))[0]
    sort(ii_1)
    points = ii_1[1:2]
    Delta_list[i] = sum(angle(Y[points]))/len(points)
    Y_Mag[i] =  abs(Y)
    Y_Angle[i] =  angle(Y)

fig = plt.figure()
ax = fig.add_subplot(111, projection = '3d')

t = linspace(-pi, pi, 1025);t = t[:-1]
fmax = 1/(t[1]-t[0])
t = t[::64]
w = linspace(-fmax*pi, fmax*pi, 65);w = w[:-1]
t, w = meshgrid(t, w)

surf = ax.plot_surface(w, t, Y_Mag.T, cmap = cm.jet, linewidth = 0, antialiased = False)
fig.colorbar(surf, shrink = 0.5, aspect = 5)
title("Magnitude Surface plot of Fragmented Chirped FT")
ylabel("Frequency")
xlabel("Time")
show()

fig = plt.figure()
ax = fig.add_subplot(111, projection = '3d')
surf = ax.plot_surface(w, t, Y_Angle.T, cmap = cm.jet, linewidth = 0, antialiased = False)
fig.colorbar(surf, shrink = 0.5, aspect = 5)
title("Phase Surface plot of Fragmented Chirped FT")
ylabel("Frequency")
xlabel("Time")
show()


#Thank you!
