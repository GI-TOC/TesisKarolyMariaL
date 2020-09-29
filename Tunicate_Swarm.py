# -*- coding: utf-8 -*-
"""
Created on Fri Sep  4 16:39:22 2020

@author: Karol
"""

import numpy as np
import math
def pobini(n, tp):
    pp = np.random.rand(tp,n)
    pp=(-5.12+(5.12+5.12)*pp)
   
    return pp

def rastrigin(pp, n, tp):
    fs=np.zeros((tp,1))
    for i in range (tp):
        sumar=0
        for j in range (n):
            sumar = sumar + ((pp[i,j]**2)-(10*math.cos(2*math.pi*(pp[i,j]))))
        sumar=sumar+10*n
        fs[i]=sumar
    
    return fs

def BEST(ppn, fsn, optimo, bestfit, tp):
    for i in range (0,tp):
        if bestfit>fsn[i]:
            optimo=ppn[i,:]
            bestfit=fsn[i]            
    return optimo, bestfit

def algebra(pp, aleatorio, tp, n, optimo):
    pp1=pp*aleatorio
    pd=np.zeros((tp,n))
    pd=pp1-optimo         
    pd=abs(pd)
    return pd

def nuevag(A, G, F, M, pmin, pmax,swarm, pp,fs):
    for i in range (Maxiteraciones):
        for j in range (2):   
            fs=rastrigin(pp, n, tp)
            c1=np.random.rand()
            c2=np.random.rand()           
            c3=np.random.rand()       
            aleatorio=np.random.rand()           
            M=M+(pmin+c1*pmax-pmin)          
            F=F+(2*c1)         
            G=G+(c2+c3)-F         
            A=A+(G/M)
            pd=algebra(pp, aleatorio, tp, n, optimo)            
            if aleatorio<=0.5:
                swarm=swarm+fs+A*pd
            else:
                swarm=swarm+fs-A*pd
       
        ppn=swarm/(2+c1)
        swarm=0
    
        i=i+1     
    for i in range (tp):
        for j in range (n):
            if ppn[i,j]<-5.12:
                ppn[i,j]=-5.12
            if ppn[i,j]>5.12:
                ppn[i,j]=5.12
        
    return ppn


#programa principal

n=int(input("ingrese el número de columnas: "))
tp=int(input("ingrese el número de agentes: "))
Maxiteraciones=int(input("ingrese el número de generaciones: "))
pp=pobini(n, tp)
fs=rastrigin(pp.copy(), n, tp)
optimo=np.zeros((1,n))+np.inf
bestfit=np.inf
optimo, bestfit=BEST(pp.copy(), fs, optimo, bestfit, tp) 
for i in range(0,Maxiteraciones):
    A=np.zeros((tp,n))
    G=np.zeros((tp,n))
    M=np.zeros((tp,n))
    F=np.zeros((tp,n))
    pmin=1
    pmax=4
    swarm=0
    fs=rastrigin(pp.copy(), n, tp)
    ppn= nuevag(A, G, F, M, pmin, pmax,swarm, pp.copy(),fs)
    fsn=rastrigin(ppn.copy(), n, tp)   
    optimo, bestfit=BEST(ppn.copy(), fsn, optimo, bestfit, tp) 
    
    pp=ppn.copy()

print(optimo, bestfit)
print("******")

    
