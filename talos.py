# -*- coding: utf-8 -*-
"""
Created on Sun Aug 26 16:22:56 2018

@author: Stanisa
"""

import numpy as np
import matplotlib.pyplot as plt
import math as math
from scipy.stats import poisson
import time

t=time.time() 

izbor= input("Choose the type of matrix!\n Type in 'cluster' or 'plankton' and press enter! ")
while izbor not in ("cluster","Cluster","CLUSTER","klaster","Plankton","PLANKTON","plankton"):
    print('Check the spelling of "cluster" or "plankton"!')
    izbor= input("Choose the type of matrix!\n Type in cluster or plankton: ")

z=int(input('Insert dimensions: '))
while z<=0:
    print('Dimensions of the grid can only be positive!')
    z=int(input('Insert dimensions: '))
    

p=int(input('Insert bacteria coverage in percent: '))
while p<=0 or p>100:
    print('Bacteria cannot cover more than 100% or less than or equal to 0% of the grid!')
    p=int(input('Insert bacteria coverage in percent: '))
    
pAntibiotik=int(input('Insert antibiotic coverage in percent: '))
while pAntibiotik<0 or pAntibiotik>100:
    print('Antibiotic cannot cover more than 100% or less than 0% of the grid!')
    pAntibiotik=int(input('Insert antibiotic coverage in percent: '))

brojGeneracija=int(input("Insert the number of generations you wish to simulate: "))


'''
dodati proveru za broj generacija da ne sme da bude manje od 0
'''


b=np.zeros((z,z))
mu=50
inClusterRatio=0.8
koncLambda=0.5
'''
podesiti da mu moze da se unese
'''
sreda =(z/2,z/2)


r_kv=((p*z*z)/math.pi/100)
rAntibiotik_kv=((pAntibiotik*z*z)/math.pi/100)

def krug(matrica,z,rSqr):
    X,Y=np.meshgrid(range(z),range(z))
    TF_mat=((X-z/2)**2+(Y-z/2)**2)<=rSqr
    P_mat=np.random.random((z,z))
    temp = (P_mat<inClusterRatio) & (TF_mat)
    b[temp]=100.0*np.random.random
    outRatio=((1-inClusterRatio)*rSqr*math.pi)/(z**2-rSqr*math.pi)
    temp2= (P_mat<outRatio) & (np.logical_not(TF_mat))
    b[temp2]=100.0*(np.random.random())
    print (b[temp2])
    
    
    
def krugAntibiotik(matrica,z,r):
    for x in range (z):
        for y in range(z):
                if ((x-z/2)**2+(y-z/2)**2)<=r :
                   matrica[x,y]=100
            
                    
def plankton(matrica,z,p):
    
     P_mat=np.random.random((z,z))
     c=(p/100)
     temp = (P_mat<c)
     matrica[temp]=100.0*np.random.random()

def generisiPlazmid(broj):
    plazmidi=[]
    for i in range (broj):
        plazmidi.append(np.random.randint(2,size=(64,)))   
    
    return plazmidi
        

def is_resistant(plasmids):
    
    res_str = np.array([0,1,0,1])
    res_begin = 4
    res_end = 8
    
    return any([all(pls[res_begin:res_end]==res_str) for pls in plasmids])


def plazmid(z):
    matrica ={'type':[],'plasmids':[],'number of plasmids':np.zeros((z,z)),'brojF':np.zeros((z,z)),'brojR':np.zeros((z,z)),'resistant':np.zeros((z,z)),'sex':np.zeros((z,z)),'lifespan':np.zeros((z,z))}
    for i in range (z):
        matrica['type'].append([])
        matrica['plasmids'].append([])
        for j in range(z):
            matrica['type'][i].append([])
            matrica['plasmids'][i].append([])
            if b[i,j]!=0:
                broj=poisson.rvs(mu)
                matrica['plasmids'][i][j] = generisiPlazmid(broj)
                matrica['number of plasmids'][i][j]=broj
                if is_resistant(matrica['plasmids'][i][j]):
                    matrica['resistant'][i][j] = 1
                for a in range (broj):
                    if ((np.random.random())<0.999):
                        matrica['type'][i][j].append('R')
                        matrica['brojR'][i][j]+=1
                    else:
                        matrica['type'][i][j].append('F')
                        matrica['brojF'][i][j]+=1
                        matrica['sex'][i][j]=1
    return matrica
        


'''
napraviti od ovih visualize funkcija
jednu fju koja prima parametar u zavisnosti od parametra
'''



def visualizePlasmids(matricaPlazmida,z):
    novaMatrica=np.zeros((z,z))
    novaMatrica=matricaPlazmida['number of plasmids']
            
    plt.style.use('classic')
    plt.matshow(novaMatrica)
    plt.title("PLASMIDS")
    plt.colorbar()
    plt.show()
        
def visualizePlasmidsR(matricaPlazmida,z):
    novaMatrica=np.zeros((z,z))
    
    novaMatrica=matricaPlazmida['brojR']
    
    plt.style.use('classic')
    plt.matshow(novaMatrica)
    plt.title("R PLASMIDS")
    plt.colorbar()
    plt.show()
    
def visualizePlasmidsF(matricaPlazmida,z):
    novaMatrica=np.zeros((z,z))
    novaMatrica=matricaPlazmida['sex']
            
    plt.style.use('classic')
    plt.matshow(novaMatrica)
    plt.title("SEX OF BACTERIA")
    plt.colorbar()
    plt.show() 
            
def visualizeBacteria(b):
    plt.style.use('classic')
    plt.matshow(b)
    plt.title("BACTERIA")
    plt.colorbar()
    plt.show()
    
def visualizeNutritives(nutrijenti):
    
    plt.style.use('grayscale')
    plt.matshow(nutrijenti)
    plt.title("NUTRITIVES")
    plt.colorbar()
    plt.show()
    
def visualizeAntibiotic(antibiotik):
    plt.style.use('classic')
    plt.matshow(antibiotik)
    plt.title("ANTIBIOTIC")
    plt.colorbar()
    plt.show()    

nbhd = np.array([[-1,-1],[-1,0],[0,-1],[1,1],[0,1],[1,0]])

def RazmnozavanjeBakterija(b,matricaPlazmida,nutrijenti):
    bTemp=np.array(b,copy=True)    
    for i in range(1,z-1):
        for  j in range(1,z-1):
          
           if(b[i,j]==0):
              # Odluci da li ces se deliti
              
              
              # Biraj gde ces se podeliti
              full_cells = []
              cells = b[i+nbhd[:,0],j+nbhd[:,1]]
              for ind,c in enumerate(cells):
                  if c:
                      full_cells.append(nbhd[ind,:])
              if full_cells:
                  ind = np.random.randint(len(full_cells))
                  bTemp[i,j]=1
                  max_food=0
                  fittest_nbhd=[]
                  fittest_nbr=[]
                  """ for n in range (len(full_cells)):
                  x_nbr=i+int((full_cells[n])[0])
                  y_nbr=j+int((full_cells[n])[1])
                  ((nutrijenti[x_nbr,y_nbr])>=max_food)
                  fittest_nbhd.append([x_nbr,y_nbr])                      
                  if n!=0:
                  fittest_nbr=fittest_nbhd[np.random.randint(0,n)]
                  matricaPlazmida['type'][i][j]=matricaPlazmida['type'][fittest_nbr[0]][fittest_nbr[1]]
                  matricaPlazmida['plasmids'][i][j]=matricaPlazmida['plasmids'][fittest_nbr[0],fittest_nbr[1]]
                  matricaPlazmida['number of plasmids'][i,j]=matricaPlazmida['number of plasmids'][fittest_nbr[0],fittest_nbr[1]]
                  matricaPlazmida['broj F'][i,j]=matricaPlazmida['broj F'][fittest_nbr[0],fittest_nbr[1]]
                  matricaPlazmida['broj R'][i,j]=matricaPlazmida['broj R'][fittest_nbr[0],fittest_nbr[1]]
                  """
                 

              
              
    b[:,:]=bTemp

if izbor=="cluster" or izbor=="Cluster" or izbor=="CLUSTER" or izbor=="klaster":
    krug(b,z,r_kv)
if izbor.upper() == "PLANKTON":
    plankton(b,z,p)
                
nutrijenti = np.zeros((z,z))
for x in range (z):
        for y in range(z):
            nutrijenti[x,y]=int(np.random.randint(0.8*brojGeneracija,1.2*brojGeneracija))

antibiotik = np.zeros((z,z))
krugAntibiotik(antibiotik,z,rAntibiotik_kv)                          
                   
matricaPlazmida=plazmid(z)


    

visualizeBacteria(b)
visualizeNutritives(nutrijenti)

visualizePlasmids(matricaPlazmida,z)

"""visualizePlasmidsR(matricaPlazmida,z)
visualizePlasmidsF(matricaPlazmida,z)
visualizeAntibiotic(antibiotik)
"""

n=brojGeneracija
while n>0:
    
    n-=1
    
    truth_vec = (b==1) & (matricaPlazmida['number of plasmids']<=45)
    nutrijenti[truth_vec]-=1
    
    truth_vec = (b==1) & (matricaPlazmida['number of plasmids']>45) & (matricaPlazmida['number of plasmids']<=55)
    nutrijenti[truth_vec]-=2
    
    truth_vec = (b==1) & (matricaPlazmida['number of plasmids']>55)
    nutrijenti[truth_vec]-=3
    truth_vec= nutrijenti<0
    nutrijenti[truth_vec]=0
    
    
    truth_vec = (b!=0) & (nutrijenti==0)
    matricaPlazmida['lifespan'][truth_vec]+=1
    
    truth_vec = (matricaPlazmida['lifespan']>3)
    b[truth_vec]=0
    matricaPlazmida['lifespan'][truth_vec]=0
    
    truth_vec=(antibiotik>60) & (matricaPlazmida['resistant']==0)
    b[truth_vec]=0
    
    RazmnozavanjeBakterija(b,matricaPlazmida,nutrijenti)
         
        
                
    visualizeNutritives(nutrijenti)
    visualizeBacteria(b)
    visualizePlasmids(matricaPlazmida,z)



"""
visualizePlasmids(matricaPlazmida,z)
visualizePlasmidsR(matricaPlazmida,z)
visualizePlasmidsF(matricaPlazmida,z)
visualizeAntibiotic(antibiotik)
"""
print('Vreme izvrsavanja ove simulacije je: ',time.time()-t)
