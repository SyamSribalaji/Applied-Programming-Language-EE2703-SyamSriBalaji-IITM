'''
EE2703-Assignment-2
Implementing Spice in Python-Part-2
Syam SriBalaji T
EE20B136
09/02/22
'''


from numpy import *
from sys import argv, exit

#1-Creating strings to identify each components
START='.circuit'
END='.end'
ACVOL_S='.ac'
AC='ac'
DC='dc'
RESISTANCE='R'
CAPACITOR='C'
INDUCTOR='L'
VOLTAGE_S='V'
CURRENT_S='I'
VCVS='E'
VCCS='G'
CCCS='F'
CCVS='H'

#2-Creating some variable which can be used in conditions
OmegaT=0
acisthere=0
dcisthere=0
acisthere1=0
dcisthere1=0
flag1=0
flag2=0

#3-Creating some empty stings to 'input and append' each lines in netlist file
NODE_A=[]
VOL_A=[]
NODE_Z=[]
VOL_Z=[]

if len(argv) != 2:#4-Checking whether the input has 2 files
    print("%s GIVE 2 VALID FILES" % argv[0])
    exit()

try:
    with open(argv[1]) as NL:#5-Opening netlist file and finding index of '.circuit','.end'&'.ac'
        lines = NL.readlines()
        start = 2; end = 1
        for line in lines:
            if START == line[:len(START)]:
                start = lines.index(line)
            elif END == line[:len(END)]:
                end = lines.index(line)
            elif ACVOL_S == line[:len(ACVOL_S)]:
            	freq = lines.index(line)
            	break
        if start >= end:#6-Condition to show error if the netlist file is invalid
            print("INVALID NETLIST FILE")
            exit(0)
        else:#7-Adding elements in the empty list which are each components(each line)
            for i in range(start+1,end):
            	if VOLTAGE_S==lines[i][:len(VOLTAGE_S)]:
            		VOL_Z=lines[i].split()
            		VOL_A.append(VOL_Z[0])
            	elif VCVS==lines[i][:len(VCVS)]:
            		VOL_Z=lines[i].split()
            		VOL_A.append(VOL_Z[0])
            	elif CCVS==lines[i][:len(CCVS)]:
            		VOL_Z=lines[i].split()
            		VOL_A.append(VOL_Z[0])
            	NODE_Z=lines[i].split()
            	NODE_A.append(NODE_Z[1])
            	NODE_A.append(NODE_Z[2])
			
except IOError:#8-Showing Error incase of invalid input netlist file
    print("INVALID NETLIST FILE")
    exit()

NODE_A=set(NODE_A)
NODE_A=(sorted)(list(NODE_A))
TOTAL_A=NODE_A+VOL_A #9-It is the list with all component names, for finding the order of nodes for each components

#10-Creating the matrices with all zero element to calculate the final result later
MATA_DC= zeros([len(TOTAL_A),len(TOTAL_A)])
MATA_AC= zeros([len(TOTAL_A),len(TOTAL_A)],complex)
MATC_DC= zeros([len(TOTAL_A),1])
MATC_AC= zeros([len(TOTAL_A),1],complex)

def findindex(k):#11-Funtion to find index of a specific component stored in list
	return(TOTAL_A.index(k))

class RES_C():#12-Class for Resistor with contructor, used to add elements of matrices of both DC and AC
	def __init__(self,res):
		res_s=res.split()
		node1=findindex(res_s[1])
		node2=findindex(res_s[2])
		value=(1/(float)(res_s[3]))
		MATA_DC[node1][node1]+=value
		MATA_DC[node1][node2]-=value
		MATA_DC[node2][node1]-=value
		MATA_DC[node2][node2]+=value
		MATA_AC[node1][node1]+=value
		MATA_AC[node1][node2]-=value
		MATA_AC[node2][node1]-=value
		MATA_AC[node2][node2]+=value

class CAP_C():#13-Class for Capacitor with contructor, used to add elements of matrices of both DC and AC
	def __init__(self,cap):
		cap_s=cap.split()
		node1=findindex(cap_s[1])
		node2=findindex(cap_s[2])
		value=1j*(OmegaT*(float)(cap_s[3]))
		MATA_AC[node1][node1]+=value
		MATA_AC[node1][node2]-=value
		MATA_AC[node2][node1]-=value
		MATA_AC[node2][node2]+=value
class IND_C():#14-Class for Inductor with contructor, used to add elements of matrices of both DC and AC
	def __init__(self,ind):
		ind_s=ind.split()
		node1=findindex(ind_s[1])
		node2=findindex(ind_s[2])
		value=1j*(-1/(float)(OmegaT*(float)(ind_s[3])))
		MATA_DC[node1][node1]+=1
		MATA_DC[node1][node2]-=1
		MATA_DC[node2][node1]-=1
		MATA_DC[node2][node2]+=1
		MATA_AC[node1][node1]+=value
		MATA_AC[node1][node2]-=value
		MATA_AC[node2][node1]-=value
		MATA_AC[node2][node2]+=value
		
class ICS_C():#15-Class for 'Independent Current Source' with contructor, used to add elements of matrix of DC
	def __init__(self,ics):
		ics_s=ics.split()
		node1=findindex(ics_s[1])
		node2=findindex(ics_s[2])
		value=(float)(ics_s[4])
		MATC_DC[node1]-=value
		MATC_DC[node2]+=value
		
class ACCS_C():#16-Class for 'AC Current Source' with contructor, used to add elements of matrix of AC
	def __init__(self,accs):
		accs_s=accs.split()
		node1=findindex(accs_s[1])
		node2=findindex(accs_s[2])
		value=(float)(accs_s[4])/2
		PHASE1=(float)(accs_s[5])
		MATC_AC[node1]-=(value*cos(pi*PHASE1/180))+1j*(value*sin(pi*PHASE1/180))
		MATC_AC[node2]+=(value*cos(pi*PHASE1/180))+1j*(value*sin(pi*PHASE1/180))
		
	
class IVS_C():#17-Class for 'Independent Voltage Source' with contructor, used to add elements of matrix of DC
	def __init__(self,ivs):
		ivs_s=ivs.split()
		node1=findindex(ivs_s[1])
		node2=findindex(ivs_s[2])
		nodez=findindex(ivs_s[0])
		value=(float)(ivs_s[4])
		MATC_DC[nodez]+=value
		MATA_DC[node1][nodez]+=1
		MATA_DC[node2][nodez]-=1
		MATA_DC[nodez][node1]+=1
		MATA_DC[nodez][node2]-=1
		
		

class ACVS_C():#18-Class for 'AC Voltage Source' with contructor, used to add elements of matrix of AC
	def __init__(self,acvs):
		acvs_s=acvs.split()
		node1=findindex(acvs_s[1])
		node2=findindex(acvs_s[2])
		nodez=findindex(acvs_s[0])
		value=(float)(acvs_s[4])/2
		PHASE2=(float)(acvs_s[5])
		MATC_AC[nodez]+=value*cos(pi*PHASE2/180)+1j*(value*sin(pi*PHASE2/180))
		MATA_AC[node1][nodez]+=1
		MATA_AC[node2][nodez]-=1
		MATA_AC[nodez][node1]+=1
		MATA_AC[nodez][node2]-=1

class VCCS_C():#19-Class for 'Voltage Controlled Current Source' with contructor, used to add elements of matrices of both AC and DC
	def __init__(self,vccs):
		vccs_s=vccs.split()
		node1=findindex(vccs_s[1])
		node2=findindex(vccs_s[2])
		node3=findindex(vccs_s[3])
		node4=findindex(vccs_s[4])
		value=(float)(vccs_s[5])
		MATA_DC[node1][node3]+=value
		MATA_DC[node2][node4]+=value
		MATA_DC[node1][node4]-=value
		MATA_DC[node2][node3]-=value
		MATA_AC[node1][node3]+=value
		MATA_AC[node2][node4]+=value
		MATA_AC[node1][node4]-=value
		MATA_AC[node2][node3]-=value
		
class VCVS_C():#20-Class for 'Voltage Controlled Voltage Source' with contructor, used to add elements of matrices of both AC and DC
	def __init__(self,vcvs):
		vcvs_s=vcvs.split()
		node1=findindex(vcvs_s[1])
		node2=findindex(vcvs_s[2])
		node3=findindex(vcvs_s[3])
		node4=findindex(vcvs_s[4])
		nodez=findindex(vcvs_s[0])
		value=(float)(vcvs_s[5])
		MATA_DC[node1][nodez]+=1
		MATA_DC[node2][nodez]-=1
		MATA_DC[nodez][node1]+=1
		MATA_DC[nodez][node2]-=1
		MATA_DC[nodez][node3]-=value
		MATA_DC[nodez][node4]+=value
		MATA_AC[node1][nodez]+=1
		MATA_AC[node2][nodez]-=1
		MATA_AC[nodez][node1]+=1
		MATA_AC[nodez][node2]-=1
		MATA_AC[nodez][node3]-=value
		MATA_AC[nodez][node4]+=value
		
class CCCS_C():#21-Class for 'Current Controlled Current Source' with contructor, used to add elements of matrices of both AC and DC
	def __init__(self,cccs):
		cccs_s=cccs.split()
		node1=findindex(cccs_s[1])
		node2=findindex(cccs_s[2])
		nodez=findindex(cccs_s[3])
		value=(float)(cccs_s[4])
		MATA_DC[node1][nodez]+=value
		MATA_DC[node2][nodez]-=value	
		MATA_AC[node1][nodez]+=value
		MATA_AC[node2][nodez]-=value	

class CCVS_C():#22-Class for 'Current Controlled Voltage Source' with contructor, used to add elements of matrices of both AC and DC
	def __init__(self,ccvs):
		ccvs_s=ccvs.split()
		node1=findindex(ccvs_s[1])
		node2=findindex(ccvs_s[2])
		nodec=findindex(ccvs_s[3])
		noded=findindex(ccvs_s[0])
		value=(float)(ccvs_s[4])
		MATA_DC[node1][noded]+=1
		MATA_DC[node2][noded]-=1
		MATA_DC[noded][node1]+=1
		MATA_DC[noded][node2]-=1
		MATA_DC[noded][nodec]-=value
		MATA_AC[node1][noded]+=1
		MATA_AC[node2][noded]-=1
		MATA_AC[noded][node1]+=1
		MATA_AC[noded][node2]-=1
		MATA_AC[noded][nodec]-=value

for i in range(start+1,end):#23-Loop to find the number of AC and DC components present in the circuit
	if VOLTAGE_S==lines[i][:len(VOLTAGE_S)]:
		k=lines[i].split()
		if k[3]=='dc':
			dcisthere1+=1
		elif k[3]=='ac':
			acisthere1+=1
	elif CURRENT_S==lines[i][:len(CURRENT_S)]:
		l=lines[i].split()
		if l[3]=='dc':
			dcisthere1+=1
		elif l[3]=='ac':
			acisthere1+=1


if acisthere1>=1:#24-Inputing the Frequency incase the circuit has AC components
	freq_s=lines[freq].split()
	FREQ=(float)(freq_s[2])
	OmegaT=(2*pi*FREQ)
			
for i in range(start+1,end):#25-Loop which goes one line by one line of the file, and creates a object for each component class(so calls the construtor and process it)
	if RESISTANCE==lines[i][:len(RESISTANCE)]:
		resistor=RES_C(lines[i])
	elif VOLTAGE_S==lines[i][:len(VOLTAGE_S)]:
		k=lines[i].split()
		if k[3]=='dc':
			indvolsour=IVS_C(lines[i])
			dcisthere+=1
		elif k[3]=='ac':
			acvolsour=ACVS_C(lines[i])
			acisthere+=1
	elif VCVS==lines[i][:len(VCVS)]:
		volconcursour=VCVS_C(lines[i])
	elif CCVS==lines[i][:len(CCVS)]:
		curconvolsour=CCVS_C(lines[i])
	elif VCCS==lines[i][:len(VCCS)]:
		volconcursour=VCCS_C(lines[i])
	elif CCCS==lines[i][:len(CCCS)]:
		curconcursour=CCCS_C(lines[i])
	elif CURRENT_S==lines[i][:len(CURRENT_S)]:
		l=lines[i].split()
		if l[3]=='dc':
			indcursour=ICS_C(lines[i])
			dcisthere+=1
		elif l[3]=='ac':
			accursour=ACCS_C(lines[i])
			acisthere+=1
	elif CAPACITOR==lines[i][:len(CAPACITOR)]:
		capacitorcall=CAP_C(lines[i])
	elif INDUCTOR==lines[i][:len(INDUCTOR)]:
		inductorcall=IND_C(lines[i])

if acisthere>=1:#26-Finding solution matrix for AC components separately
	ind=findindex('GND')
	MATA_AC[ind][ind]+=1
	MATB_AC=linalg.solve(MATA_AC,MATC_AC)
	flag1+=1

if dcisthere>=1:#27-Finding solution matrix for DC components separately
	ind=findindex('GND')
	MATA_DC[ind][ind]+=1
	MATB_DC=linalg.solve(MATA_DC,MATC_DC)
	flag2+=1

#28-Loops to display the final output(i.e. The Voltage in each node and Current thourgh voltage sources)	
if (flag1==1)&(flag2==1):
	for i in range(0,len(NODE_A)):
		print("The voltage of the node %d is %s + Re(%s exp(j* %f *t))" %(i+1,MATB_DC[i],MATB_AC[i],OmegaT))
	for i in range(len(NODE_A),len(TOTAL_A)):
		print("The current through the voltage source %s is %s + Re(%s exp(j* %f *t))" %(TOTAL_A[i],MATB_DC[i],MATB_AC[i],OmegaT))
elif (flag1==1)&(flag2!=1):
	for i in range(0,len(NODE_A)):
		print("The voltage of the node n%d is Re(%s exp(j* %f *t))" %(i+1,MATB_AC[i],OmegaT))
	for i in range(len(NODE_A),len(TOTAL_A)):
		print("The current through the voltage source %s is Re(%s exp(j* %f *t))" %(TOTAL_A[i],MATB_AC[i],OmegaT))
elif (flag1!=1)&(flag2==1):
	for i in range(0,len(NODE_A)):
		print("The voltage of the node n%d is %s" %(i+1,MATB_DC[i]))
	for i in range(len(NODE_A),len(TOTAL_A)):
		print("The current through the voltage source %s is %s" %(TOTAL_A[i],MATB_DC[i]))

	
	
	
	
		
		
		

