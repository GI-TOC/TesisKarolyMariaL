# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 10:12:02 2020

@author: Karol
"""
import numpy as np


#pobini
def pobini(tp,n):
    matriz=np.zeros((tp,n))
    for i in range (tp):
        vector=np.random.permutation(n)+1
        matriz[i,:]=vector
    return matriz

#evaluación
def evaluacion(matriz,tp,n,dem,cap):
   ruta=np.zeros((tp,(n*2)+1))
   cont=0
   #asignacion primer cliente
   for i in range (tp):
       ruta[i,1]=matriz[i,0]
       #asignacion del segundo cliente hasta el ultimo
   for i in range (tp):
       cont=dem[matriz[i,0]-1]
       pos=1
       j=1
       while j<n:           
           cont = cont + dem[matriz[i,j]-1]
           if cont <= cap:
               ruta[i,pos+1] = matriz[i,j]
               pos=pos+1
               j=j+1               
           else:
               ruta[i,pos+1]=0
               pos=pos+1
               cont=0
   return ruta

#costo_matriz
def costo(a,n,tp,ruta2,dem):
    arco=np.array([[0,1,2,1,1,2],[1,0,1,1,2,2],[1,1,0,1,2,1],[1,2,2,0,1,1],[2,1,1,2,0,2],[1,2,1,1,2,0]])
    m=(n*2)-1
    valor=np.zeros((tp,1))   
    cont1=0
    for i in range (tp):
        cont1=0
        for j in range (m):
            cont1 = cont1 + arco[ruta2[i,j],ruta2[i,j+1]]
            valor[i,0]=cont1    
    return valor
#torneo
def torneo(n, tp, a, valor):
    padres=np.zeros((tp,n))
    for i in range (tp):
        ale1=np.random.randint(0,tp-1)
        ale2=np.random.randint(0,tp-1)
        while ale1==ale2:
            ale1=np.random.randint(0,tp-1)
            ale2=np.random.randint(0,tp-1)
        if valor[ale1]<= valor[ale2]:
            padres[i,:]=a[ale1,:]
        else:
            padres[i,:]=a[ale2,:]    
    return padres
#SB2OX
def cruzamiento(n,tp,padres, pc):
    hijos=np.zeros((tp,n))
    for i in range (0,tp,2):
        al1=np.random.rand()
        if al1<=pc:
            al2=np.random.randint(0,n-1)
            al3=np.random.randint(0,n-1)            
            while al2>=al3:
                al2=np.random.randint(0,n-1)
                al3=np.random.randint(0,n-1)
            #print("ale2")
            #print(al2)
            #print("ale3")
            #print(al3)
            for j in range (n-1):
                if padres[i,j]==padres[i+1,j] and padres[i,j+1]==padres[i+1,j+1]:
                    hijos[i,j]=padres[i,j]
                    hijos[i+1,j]=padres[i+1,j]
                    hijos[i,j+1]=padres[i,j+1]
                    hijos[i+1,j+1]=padres[i+1,j+1]
                    
            for j in range (n):                
                a=padres[i,al2:al3+1]
                c=padres[i+1,al2:al3+1]
                hijos[i,al2:al3+1]=a
                hijos[i+1,al2:al3+1]=c              
            
                 
        else:
            hijos[i,:]=padres[i,:]
            hijos[i+1,:]=padres[i+1,:]           
            
    return hijos
def cruce2(hijos,padres,tp,n):
    for i in range (0,tp,2):
        temp=0
        cont=0
        temp1=0
        cont1=0
        for j in range (n):
            temp = np.where(hijos[i]==padres[i+1,j])            
            if temp[0].size<=0:
                cont=np.where(hijos[i]==0)            
                hijos[i, cont[0][0]]=padres[i+1,j]             

            temp1 = np.where(hijos[i+1]==padres[i,j])    
            if temp1[0].size<=0:             
                cont1=np.where(hijos[i+1]==0)
                hijos[i+1,cont1[0][0]]=padres[i,j]
    return hijos
def mutacion(hijos,n,tp,pm):
    hijosm=hijos.copy()
    a=0
    b=0
    for i in range (tp):
        all1=np.random.rand()
        if all1<=pm:
            ale1=np.random.randint(0,n)
            ale2=np.random.randint(0,n)
            while ale1==ale2:
                ale1=np.random.randint(0,n)
                ale2=np.random.randint(0,n)
            a=hijosm[i,ale1]
            b=hijosm[i,ale2]
            hijosm[i,ale1]=b
            hijosm[i,ale2]=a        
    return hijosm
def actualizar(hijos, valorh):
    optimo=matriz[0,:]
    optfit=valor[0]
    for i in range (0,tp):
        if optfit>valorh[i]:
            optimo=hijos[i,:]
            optfit=valorh[i]
    
    return optimo, optfit

    



#Programaprincipal
tp=int(input("ingrese el tamaño de la población: "))
n=5
if tp%2==0:
    tp=tp
else:
    tp=tp+1
pc=float(input("ingrese la probabilidad de cruzamiento: "))
pm=float(input("ingrese la probabilidad de mutación: "))
ngen=int(input("ingrese el número de generaciones: "))
matriz=pobini(tp,n)
matriz=matriz.astype(int)
#print(matriz)
dem=[5, 7, 10, 9, 11]
cap=25  
for i in range(0,ngen):
    ruta=evaluacion(matriz,tp,n,dem,cap)
    ruta=ruta.astype(int)
    #print("la rutas de la población inicial")
    #print(ruta)
    valor=costo(matriz,n,tp,ruta,dem)
    valor=valor.astype(int)
    #print("El costo de la ruta es: ")
    #print(valor)
    padres=torneo(n, tp, matriz, valor)
    padres=padres.astype(int)
    #print("Matriz de padres")
    #print(padres)
    hijos=cruzamiento(n,tp,padres,pc)
    hijos=hijos.astype(int)
    #print("matriz de hijos")
    #print(hijos)
    hijos=cruce2(hijos,padres,tp,n)
    #print("******")
    #print(hijos)
    rutah=evaluacion(hijos,tp,n,dem,cap)
    valorh=costo(hijos,n,tp,ruta,dem)
    valorh=valorh.astype(int)
    optimo, optfit=actualizar(hijos, valorh)
    matriz=hijos.copy()
print(optimo)
print(optfit)




