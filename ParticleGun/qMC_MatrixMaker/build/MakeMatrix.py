import numpy as np

code=input("Input Code name for material:")
bins=np.arange(0,20001,5)
#---------------------
#gXg

Matrix=[]
with open("./RunTest.txt","r") as fi:
    id = []
    for ln in fi:
        if ln.startswith("G:"):
            id.append(ln[2:])
L=len(id)
######            
NperE=1000
######



data=[]
for i in range(L//(2*NperE),L//NperE):
    data.append(id[NperE*i:NperE*(i+1)])
    
for energy in data:
    eOut=[]    
    for event in energy:
        event=event.split(',')
        event=[i for i in event if i!='\n']        
        event=[float(i)*1e3 for i in event if i!='0']    
        eOut+=event
    
    eOut=np.array(eOut)
    if sum(np.histogram(eOut, bins)[0])==0:
        s_norm=1
    else:
        s_norm=sum(np.histogram(eOut, bins)[0])
        
    prob=(eOut.size/NperE) * np.histogram(eOut, bins)[0] / s_norm
    #print(out_spect)    
    Matrix.append(prob)

fileMat=open("./matrices/Results_g"+code+"g.txt","w")

for i in Matrix:
    for j in i:
        fileMat.write(str(j)+"\t")        
    fileMat.write("\n")

fileMat.close()

#---------------------
#gXn

Matrix=[]
with open("./RunTest.txt","r") as fi:
    id = []
    for ln in fi:
        if ln.startswith("N:"):
            id.append(ln[2:])
L=len(id)
######            
NperE=1000
######



data=[]
for i in range(L//(2*NperE),L//NperE):
    data.append(id[NperE*i:NperE*(i+1)])
    
for energy in data:
    eOut=[]    
    for event in energy:
        event=event.split(',')
        event=[i for i in event if i!='\n']        
        event=[float(i)*1e3 for i in event if i!='0']    
        eOut+=event
    
    eOut=np.array(eOut)
    if sum(np.histogram(eOut, bins)[0])==0:
        s_norm=1
    else:
        s_norm=sum(np.histogram(eOut, bins)[0])
        
    prob=(eOut.size/NperE) * np.histogram(eOut, bins)[0] / s_norm
    #print(out_spect)    
    Matrix.append(prob)

fileMat=open("./matrices/Results_g"+code+"n.txt","w")

for i in Matrix:
    for j in i:
        fileMat.write(str(j)+"\t")        
    fileMat.write("\n")

fileMat.close()


#---------------------
#nXg

Matrix=[]
with open("./RunTest.txt","r") as fi:
    id = []
    for ln in fi:
        if ln.startswith("G:"):
            id.append(ln[2:])
L=len(id)
######            
NperE=1000
######



data=[]
for i in range(L//(2*NperE)):
    data.append(id[NperE*i:NperE*(i+1)])
    
for energy in data:
    eOut=[]    
    for event in energy:
        event=event.split(',')
        event=[i for i in event if i!='\n']        
        event=[float(i)*1e3 for i in event if i!='0']    
        eOut+=event
    
    eOut=np.array(eOut)
    if sum(np.histogram(eOut, bins)[0])==0:
        s_norm=1
    else:
        s_norm=sum(np.histogram(eOut, bins)[0])
        
    prob=(eOut.size/NperE) * np.histogram(eOut, bins)[0] / s_norm
    #print(out_spect)    
    Matrix.append(prob)

fileMat=open("./matrices/Results_n"+code+"g.txt","w")

for i in Matrix:
    for j in i:
        fileMat.write(str(j)+"\t")        
    fileMat.write("\n")

fileMat.close()

#---------------------
#nXn

Matrix=[]
with open("./RunTest.txt","r") as fi:
    id = []
    for ln in fi:
        if ln.startswith("N:"):
            id.append(ln[2:])
L=len(id)
######            
NperE=1000
######



data=[]
for i in range(L//(2*NperE)):
    data.append(id[NperE*i:NperE*(i+1)])
    
for energy in data:
    eOut=[]    
    for event in energy:
        event=event.split(',')
        event=[i for i in event if i!='\n']        
        event=[float(i)*1e3 for i in event if i!='0']    
        eOut+=event
    
    eOut=np.array(eOut)
    if sum(np.histogram(eOut, bins)[0])==0:
        s_norm=1
    else:
        s_norm=sum(np.histogram(eOut, bins)[0])
        
    prob=(eOut.size/NperE) * np.histogram(eOut, bins)[0] / s_norm
    #print(out_spect)    
    Matrix.append(prob)

fileMat=open("./matrices/Results_n"+code+"n.txt","w")

for i in Matrix:
    for j in i:
        fileMat.write(str(j)+"\t")        
    fileMat.write("\n")

fileMat.close()

