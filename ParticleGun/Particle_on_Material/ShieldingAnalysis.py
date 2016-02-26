import numpy as np

#Build Matrices from .txt files
print("Building Matrices...")
#Open Pb Matrices
#nPbn
file_nPbn=open('./Results_nPbn.txt')
M_nPbn=[]
for line in file_nPbn:
    line.strip()
    row=[float(elt) for elt in line.split()]
    M_nPbn.append(row)
M_nPbn=np.array(M_nPbn)
print ("nPbn Ready!")

#nPbg
file_nPbg=open('./Results_nPbg.txt')
M_nPbg=[]
for line in file_nPbg:
    line.strip()
    row=[float(elt) for elt in line.split()]
    M_nPbg.append(row)
M_nPbg=np.array(M_nPbg)
print ("nPbg Ready!")

#gPbn
file_gPbn=open('./Results_gPbn.txt')
M_gPbn=[]
for line in file_gPbn:
    line.strip()
    row=[float(elt) for elt in line.split()]
    M_gPbn.append(row)
M_gPbn=np.array(M_gPbn)
print ("gPbn Ready!")

#gPbg
file_gPbg=open('./Results_gPbg.txt')
M_gPbg=[]
for line in file_gPbg:
    line.strip()
    row=[float(elt) for elt in line.split()]
    M_gPbg.append(row)
M_gPbg=np.array(M_gPbg)
print ("gPbg Ready!")

#Create Pb Material Matrix
M_Pb=np.array([[M_nPbn,M_nPbg],[M_gPbn,M_gPbg]])

print("Pb Matrix Complete.")

#Open Concrete Matrices

#nConn
file_nConn=open('./Results_nConn.txt')
M_nConn=[]
for line in file_nConn:
    line.strip()
    row=[float(elt) for elt in line.split()]
    M_nConn.append(row)
M_nConn=np.array(M_nConn)
print ("nConn Ready!")

#nCong
file_nCong=open('./Results_nCong.txt')
M_nCong=[]
for line in file_nCong:
    line.strip()
    row=[float(elt) for elt in line.split()]
    M_nCong.append(row)
M_nCong=np.array(M_nCong)
print ("nCong Ready!")

#gConn
file_gConn=open('./Results_gConn.txt')
M_gConn=[]
for line in file_gConn:
    line.strip()
    row=[float(elt) for elt in line.split()]
    M_gConn.append(row)
M_gConn=np.array(M_gConn)
print ("gConn Ready!")

#gCong
file_gCong=open('./Results_gCong.txt')
M_gCong=[]
for line in file_gCong:
    line.strip()
    row=[float(elt) for elt in line.split()]
    M_gCong.append(row)
M_gCong=np.array(M_gCong)
print ("gCong Ready!")

#Create Water Material Matrix
M_Con=np.array([[M_nConn,M_nCong],[M_gConn,M_gCong]])
print("Concrete Matrix Complete.")
#Open Water Matrices

#nWatn
file_nWatn=open('./Results_nWatn.txt')
M_nWatn=[]
for line in file_nWatn:
    line.strip()
    row=[float(elt) for elt in line.split()]
    M_nWatn.append(row)
M_nWatn=np.array(M_nWatn)
print ("nWatn Ready!")

#nWatn
file_nWatg=open('./Results_nWatg.txt')
M_nWatg=[]
for line in file_nWatg:
    line.strip()
    row=[float(elt) for elt in line.split()]
    M_nWatg.append(row)
M_nWatg=np.array(M_nWatg)
print ("nWatg Ready!")

#gWatn
file_gWatn=open('./Results_gWatn.txt')
M_gWatn=[]
for line in file_gWatn:
    line.strip()
    row=[float(elt) for elt in line.split()]
    M_gWatn.append(row)
M_gWatn=np.array(M_gWatn)
print ("gWatn Ready!")

#gWatg
file_gWatg=open('./Results_gWatg.txt')
M_gWatg=[]
for line in file_gWatg:
    line.strip()
    row=[float(elt) for elt in line.split()]
    M_gWatg.append(row)
M_gWatg=np.array(M_gWatg)
print ("gWatg Ready!")

#Create Water Material Matrix
M_Wat=np.array([[M_nWatn,M_nWatg],[M_gWatn,M_gWatg]])
print("Water Matrix Complete.")

#Open Polyethylene Matrices

#nPoln
file_nPoln=open('./Results_nPoln.txt')
M_nPoln=[]
for line in file_nPoln:
    line.strip()
    row=[float(elt) for elt in line.split()]
    M_nPoln.append(row)
M_nPoln=np.array(M_nPoln)
print ("nPoln Ready!")

#nPolg
file_nPolg=open('./Results_nPolg.txt')
M_nPolg=[]
for line in file_nPolg:
    line.strip()
    row=[float(elt) for elt in line.split()]
    M_nPolg.append(row)
M_nPolg=np.array(M_nPolg)
print ("nPolg Ready!")

#gPoln
file_gPoln=open('./Results_gPoln.txt')
M_gPoln=[]
for line in file_gPoln:
    line.strip()
    row=[float(elt) for elt in line.split()]
    M_gPoln.append(row)
M_gPoln=np.array(M_gPoln)
print ("gPoln Ready!")

#gPolg
file_gPolg=open('./Results_gPolg.txt')
M_gPolg=[]
for line in file_gPolg:
    line.strip()
    row=[float(elt) for elt in line.split()]
    M_gPolg.append(row)
M_gPolg=np.array(M_gPolg)
print ("gPolg Ready!")

#Create Polyethylene Material Matrix
M_Pol=np.array([[M_nPoln,M_nPolg],[M_gPoln,M_gPolg]])
print("Polyethylene Matrix Complete.")

#nGdn
file_nGdn=open('./Results_nGdn.txt')
M_nGdn=[]
for line in file_nGdn:
    line.strip()
    row=[float(elt) for elt in line.split()]
    M_nGdn.append(row)
M_nGdn=np.array(M_nGdn)
print ("nGdn Ready!")

#nGdg
file_nGdg=open('./Results_nGdg.txt')
M_nGdg=[]
for line in file_nGdg:
    line.strip()
    row=[float(elt) for elt in line.split()]
    M_nGdg.append(row)
M_nGdg=np.array(M_nGdg)
print ("nGdg Ready!")

#gGdn
file_gGdn=open('./Results_gGdn.txt')
M_gGdn=[]
for line in file_gGdn:
    line.strip()
    row=[float(elt) for elt in line.split()]
    M_gGdn.append(row)
M_gGdn=np.array(M_gGdn)
print ("gGdn Ready!")

#gGdg
file_gGdg=open('./Results_gGdg.txt')
M_gGdg=[]
for line in file_gGdg:
    line.strip()
    row=[float(elt) for elt in line.split()]
    M_gGdg.append(row)
M_gGdg=np.array(M_gGdg)
print ("gGdg Ready!")

#Create Gadolinium Material Matrix
M_Gd=np.array([[M_nGdn,M_nGdg],[M_gGdn,M_gGdg]])
print("Gadolinium Matrix Complete.")

#nGen
file_nGen=open('./Results_nGen.txt')
M_nGen=[]
for line in file_nGen:
    line.strip()
    row=[float(elt) for elt in line.split()]
    M_nGen.append(row)
M_nGen=np.array(M_nGen)
print ("nGen Ready!")

#nGeg
file_nGeg=open('./Results_nGeg.txt')
M_nGeg=[]
for line in file_nGeg:
    line.strip()
    row=[float(elt) for elt in line.split()]
    M_nGeg.append(row)
M_nGeg=np.array(M_nGeg)
print ("nGeg Ready!")

#gGen
file_gGen=open('./Results_gGen.txt')
M_gGen=[]
for line in file_gGen:
    line.strip()
    row=[float(elt) for elt in line.split()]
    M_gGen.append(row)
M_gGen=np.array(M_gGen)
print ("gGen Ready!")

#gGeg
file_gGeg=open('./Results_gGeg.txt')
M_gGeg=[]
for line in file_gGeg:
    line.strip()
    row=[float(elt) for elt in line.split()]
    M_gGeg.append(row)
M_gGeg=np.array(M_gGeg)
print ("gGeg Ready!")

#Create Germanium Material Matrix
M_Ge=np.array([[M_nGen,M_nGeg],[M_gGen,M_gGeg]])
print("Germanium Matrix Complete.")

#nSin
file_nSin=open('./Results_nSin.txt')
M_nSin=[]
for line in file_nSin:
    line.strip()
    row=[float(elt) for elt in line.split()]
    M_nSin.append(row)
M_nSin=np.array(M_nSin)
print ("nSin Ready!")

#nSig
file_nSig=open('./Results_nSig.txt')
M_nSig=[]
for line in file_nSig:
    line.strip()
    row=[float(elt) for elt in line.split()]
    M_nSig.append(row)
M_nSig=np.array(M_nSig)
print ("nSig Ready!")

#gSin
file_gSin=open('./Results_gSin.txt')
M_gSin=[]
for line in file_gSin:
    line.strip()
    row=[float(elt) for elt in line.split()]
    M_gSin.append(row)
M_gSin=np.array(M_gSin)
print ("gSin Ready!")

#gSig
file_gSig=open('./Results_gSig.txt')
M_gSig=[]
for line in file_gSig:
    line.strip()
    row=[float(elt) for elt in line.split()]
    M_gSig.append(row)
M_gSig=np.array(M_gSig)
print ("gSig Ready!")

#Create Germanium Material Matrix
M_Si=np.array([[M_nSin,M_nSig],[M_gSin,M_gSig]])
print("Silicon Matrix Complete.")

print(" All Matrices Built.")



#Input Spectrum of Neutrons
from scipy.stats import norm
E = np.arange(10,10001,10)
#I_spec_n=1000000000*maxwell.pdf((E+2000)/1.4,0,2000)
#Gaussian Spectrum
I_spec_n=(10**6)*norm.pdf(E,100,100)
#I_spec_n=np.array(1000*[0])
print("Input Spectrum of Neutrons Built")


#Input Spectrum of Gammas
#from scipy.stats import norm
E = np.arange(10,10001,10)
#Flat Spectrum
I_spec_g=np.array(1000*[2000])
#I_spec_g=(100000000/E)+(10**3)
print("Input Spectrum of Gamma Built")

#Input Spectrum
#I_spec=np.array([I_spec_n,I_spec_g])

#Compute Output Spectrum based length in cm
def Layers(I_n,I_g,M,l):    
    O_n,O_g=I_n,I_g
    for i in range(l):
        O_n , O_g = M[0][0].transpose().dot(O_n)+M[1][0].transpose().dot(O_g) , M[0][1].transpose().dot(O_n)+M[1][1].transpose().dot(O_g)
    return O_n,O_g
 

#Compute Output Spectra for Shielding Layers based on Input Scheme
print("Build the Shielding scheme you wish to analyse")

typs=[]
thicks=[]
flag=True
dict={'1':"Pb",'2':"HDC",'3':"Water",'4':"Poly",'5':"Gd","6":"Ge","7":"Si"}
while(flag):
  typ=input("Add a Layer\n Press 1 for Pb \n Press 2 for HDC \n Press 3 for Water \n Press 4 for Polyethylene \n Press 5 for Gadolinium \n Press 6 for Germanium \n Press 7 for Silicon \nInput: ")
  thick=input("Input thickness in cm: ")
  typs.append(dict[typ])
  thicks.append(thick)
  
  inp=input("Press 0 if done, Press Enter/Return to add another Layer\n")
  if inp=="0":
      flag=False


O_spec_n,O_spec_g=I_spec_n,I_spec_g
for i in range(len(typs)):
    mater={'Pb':M_Pb,'HDC':M_Con,'Water':M_Wat, 'Poly': M_Pol,'Gd': M_Gd,'Ge':M_Ge,'Si':M_Si}
    O_spec_n,O_spec_g=Layers(O_spec_n,O_spec_g,(mater[typs[i]]),int(thicks[i]))

Details="-->|"
for i in range(len(typs)):
    Details+=typs[i]+" ("+thicks[i]+" cm)|"

print("Shielding Scheme:",Details)
import matplotlib.pyplot as plt
plt.figure(figsize=(12,12))

#Plot Input Spectrum of Neutrons
i_n,=plt.plot(E,I_spec_n, color='b')

#Plot Input Spectrum of Gammas
import matplotlib.pyplot as plt
i_g,=plt.plot(E,I_spec_g, color='y')

#Plot Output Spectrum of Neutrons
import matplotlib.pyplot as plt
o_n,=plt.plot(E,O_spec_n, color='b', linestyle="--")

#Plot Output Spectrum of Gammas
import matplotlib.pyplot as plt
o_g,=plt.plot(E,O_spec_g, color='y', linestyle="--")


plt.xscale("log")
plt.title("Shielding\n"+Details)
plt.xlabel("Energy (keV)")
plt.ylabel("Number of Entries")
plt.legend([i_n,i_g,o_n,o_g], ["Input Neutron Spectrum","Input Gamma Spectrum","Output Neutron Spectrum","Output Gamma Spectrum"])
#plt.yscale("log")
