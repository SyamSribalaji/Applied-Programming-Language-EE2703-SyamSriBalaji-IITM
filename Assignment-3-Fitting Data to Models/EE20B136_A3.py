'''
EE2703-Assignment-3
Fitting Data to Models
Syam SriBalaji T
EE20B136
18/02/22
'''

from pylab import*
from scipy.special import*
from numpy import*
exec(open("generate_data.py").read())
show()
Data=loadtxt('fitting.dat')
Time=Data[:,0]
A=1.05
B=-0.105

#***********************************************************************************************
#Question-4-(i)
#Here we find and plot the correct values of the graph needed(i.e. without any noise)

def Truedata():
	x=Time
	A=1.05
	B=-0.105
	y=A*jn(2,x)+B*x
	plot(x,y)
	xlabel('Time')
	ylabel('f(t)=y=1.05*jn(2,x)+(-0.105)*x')
	title('Q4: Correct data to be fitted in graph')
	show()
	return y
TV=zeros([len(Time),1])
T1=Truedata()
for d in range(0,len(Time)):
	TV[d][0]+=T1[d]

#***********************************************************************************************
#Question-4-(ii)
#Here we input and plot the each column's data with different Standard deviation

for u in range(0,k):
	y=Data[:,u+1]
	x=Time
	plot(x,y,label='σ%d=%f'%(u+1,scl[u]))
x=Time
y=A*jn(2,x)+B*x
plot(x,y,label='True values')
xlabel('Time')
ylabel('Datas in each column')
title('Q4: Noised data graph for each Standard deviation')
legend()
show()

#***********************************************************************************************
#Question-5
#Here we plot the Error bars for the datas in column-1(with σ=0.10)

Data1=Data[:,1]
stdev=scl[0]
x=Time
y=A*jn(2,x)+B*x
errorbar(Time[::5],Data1[::5],stdev,fmt='ro',label='Error bar for Column-1 datas')
plot(x,y,label='True value graph')
xlabel('Time')
ylabel('Column-1-datas')
title('Q5: Data points for σ=0.10(in Error bars) along with exact function')
legend()
show()

#***********************************************************************************************
#Question-6
#Here we create the True value graph in different method and verify it with the answer got in Question-4

def Matrix_Truevalue(A,B):
	M=zeros([len(Time),2])
	P=zeros([2,1])
	x=Time
	P[0][0]+=A
	P[1][0]+=B
	M[:,0]+=jn(2,x)
	M[:,1]+=x
	G=dot(M,P)
	print('\n\nCorrect Datas which need to be plotted')
	print(G)
	return(G)
G1=Matrix_Truevalue(1.05,-0.105)
Dif_M=zeros([len(Time),1])
Dif_M=G1-TV
D=std(Dif_M[:,0])
print('\n\nThe value of D is = %f'%D)
E=(int)(D)
if(E==0):
	print('Both the Matrices are Same')
else:
	print('Both the Matrices are not same')

#***********************************************************************************************
#Question-7
#Here we find Mean squared error matrix by taking different A and B values

N=21
Avalues=linspace(0,2,N)
Bvalues=linspace(-0.2,0,N)
Err_val= zeros([N,N])
def Er_Matrix1():
	for i in range(0,N):
		for j in range(0,N):
			G_value=Avalues[i]*jn(2,Time)+Bvalues[j]*Time
			for k in range(0,len(Time)):
				Err_val[i][j]+=(1/101)*((Data1[k]-G_value[k])**2)
	print('\n\nMean squared error matrix for Datas of Column1')	
	print(Err_val)
	return(Err_val)
Z_Er=Er_Matrix1()

#***********************************************************************************************
#Question-8
#Here we plot the contour graph using the Mean squared matrix found in Question-7

xlabel('A')
ylabel('B')
title('Q8: Contour plot of ε ij')
CS=contour(Avalues,Bvalues,Z_Er,20)
plot(A, B, marker="o", markersize=10,label='Exact value')
clabel(CS,inline=1, fontsize=10)
legend()
show()

#***********************************************************************************************
#Question-9
#Here we find best estimated values of A and B for the column-1(with σ=0.10)

M=zeros([len(Time),2])
M[:,0]+=jn(2,Time)
M[:,1]+=Time
y=Data1
p, res, rnk, s = lstsq(M, y,rcond=None)
print('\n\nThe best estimated values of A and B for Datas of Column1')
print(p)
print('\n\n')

#***********************************************************************************************
#Question-10
#Here we run 'generate_data.py' multiple times and create different sets of noised datas. And find and plot errrors in estimation of A and B
 
NOI=10
A1=zeros([NOI,k,2])
for i in range(0,NOI):
	exec(open("generate_data.py").read())
	Datax=loadtxt('fitting.dat')
	Time=Datax[:,0]
	M=zeros([len(Time),2])
	M[:,0]+=jn(2,Time)
	M[:,1]+=Time
	p= zeros([k,2])
	for s in range(0,k):
		y=Datax[:,s+1]
		A1[i,s,:], res, rnk, s = lstsq(M, y,rcond=None)
show()

A_Err=list(0 for c in range(0,k))
B_Err=list(0 for c in range(0,k)) 
for j in range(0,k):
	A_Err[j]=std(A1[:,j,0])
	B_Err[j]=std(A1[:,j,1])
x=scl
y1=A_Err
y2=B_Err
plot(x,y1,marker = 'o',linestyle='dashed',color='red',label='Error in A')
plot(x,y2,marker = 'o',linestyle='dashed',color='green',label='Error in B')
xlabel('Noise standard deviation')
ylabel('MS error')
title('Q10: Variation of error with noise in A & B')
legend()
show()

#***********************************************************************************************
#Question-11
#Here we do the same process from Question-10 but now by using loglog()function

loglog(x,y1,'ro',label='Error in A')
stem(x,y1,'ro--')
loglog(x,y2,'go',label='Error in B')
stem(x,y2,'go--')
xlabel('σn')
ylabel('MS error')
title('Q11: Variation of error with noise in A & B (by using loglog())')
legend()
show()

		
#PLEASE NOTE: I have commented lines used for plotting and displaying graph in 'generate_data.py' file. As it is done is this file itself. Thank you!

