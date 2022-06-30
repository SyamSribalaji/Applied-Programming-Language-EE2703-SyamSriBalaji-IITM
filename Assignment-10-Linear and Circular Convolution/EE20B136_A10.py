'''
EE2703-Assignment-10
Linear and Circular Convolution
Syam SriBalaji T
EE20B136
11/05/22
'''

from numpy import *
from pylab import *
from scipy import signal
import argparse 
import csv
import time


#**********************************************************************************************************
#Question-1
#Reading and Receiving datas from the given file

parser = argparse.ArgumentParser()
parser.add_argument('--file1',type=str,default='./h.csv')
parser.add_argument('--file2',type=str,default='./x1.csv')
args = parser.parse_args()

h = genfromtxt(args.file1,delimiter=',')
w, H = signal.freqz(h)

#**********************************************************************************************************
#Question-2
#Magnitude and Phase plot of the filter

subplot(2,1,1)
plot(w,20*log10(abs(H)),'g')
ylabel(r'$|Y|$')
title("Magnitude and Phase response for the given Filter")
subplot(2,1,2)
plot(w,angle(H),'g')
ylabel("$∠Y$")
xlabel("$\omega$")
show()

#**********************************************************************************************************
#Question-3
#Graphing for x = cos(0.2*π*n) + cos(0.85*π*n)

n = linspace(1,2**10,2**10)
x = cos(0.2*pi*n) + cos(0.85*pi*n)
plot(n,x,'r')
xlabel("n")
ylabel("x")
title("Plot of $x = \cos(0.2\pi t)+\cos(0.85\pi t)$")
xlim([1,100])
show()

#**********************************************************************************************************
#Question-4
#Linear convolution of x = cos(0.2*π*n) + cos(0.85*π*n)

y =  convolve(x,h)
plot(range(len(n)+len(h)-1),y,'r')
xlabel("n")
ylabel("y")
title("Plot of $y=x ⊗ h$ using Linear convolution")
xlim([1,100])
show()

#**********************************************************************************************************
#Question-5
#Circular convolution of y = cos(0.2*π*n) + cos(0.85*π*n) using DFTs techniques

x_ = concatenate((x,zeros(len(h)-1)))
y1= np.fft.ifft(np.fft.fft(x_) * np.fft.fft( concatenate( (h,zeros(len(x_)-len(h))) )))
plot(range(len(y1)),y1,'r')
xlabel("n")
ylabel("y")
title("Plot of $y=x ⊗ h$ using DFTs techniques")
xlim([1,100])
show()

#**********************************************************************************************************
#Question-6
#Linear convolution using Circular convolution for x = cos(0.2*π*n) + cos(0.85*π*n)

def Circ_conv(x,h):
    P1 = len(h)
    n_ = int(ceil(log2(P1)))
    h_ = concatenate((h,zeros(int(2**n_)-P1)))
    P2 = len(h_)
    n1 = int(ceil(len(x)/2**n_))
    x_ = concatenate((x,zeros(n1*(int(2**n_))-len(x))))
    y = zeros(len(x_)+len(h_)-1)
    for i in range(n1):
        temp = concatenate((x_[i*P2:(i+1)*P2],zeros(P2-1)))
        y[i*P2:(i+1)*P2+P2-1] += np.fft.ifft(np.fft.fft(temp) * np.fft.fft( concatenate( (h_,zeros(len(temp)-len(h_))) ))).real
    return y

Y = Circ_conv(x,h)
plot(range(len(Y)),Y,'r')
xlabel("n")
ylabel("y")
title("Plot of $y=x ⊗ h$ using Circular convolution")
xlim([0,100])
show()

#**********************************************************************************************************
#Question-7
#Correlation output plot of Zadoff-Chu sequences

lines = []
with open(args.file2,'r') as file2:
    csvreader = csv.reader(file2)
    for row in csvreader:
        lines.append(row)
lines2 = []
for line in lines:
    line = list(line[0])
    try :
        line[line.index('i')]='j'
        lines2.append(line)
    except ValueError:
        lines2.append(line)
        continue

x = [complex(''.join(line)) for line in lines2]
X = np.fft.fft(x)
x2 = roll(x,5)
noc = np.fft.ifftshift(correlate(x2,x,'full'))
print("The number of correlation is %s" %(len(noc)))

stem(linspace(0,len(noc)-1,len(noc)),abs(noc),'g','go','r')
xlabel("t")
ylabel("Correlation")
title("Correlation plot of Zadoff-Chu sequences")
xlim([1,20])
show()

#Thank you!



