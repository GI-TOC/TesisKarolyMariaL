# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 22:56:22 2019

@author: Helman Hernández Riaño
"""
# ALGORITMO GENÉTICO SIMPLE PARA LA MINIMIZACIÓN DE LA
# FUNCIÓN RASTRIGIN

import numpy as np
#import os as osi
#osi.system("clear")

# Función de Población Inicial
def PobIni(tp,n,a,b):
    PI=a+(b-a)*np.random.rand(tp,n)
    return PI

# Función de evaluación (fitness)
def fitness(PI,tp,n,A):
    fit=np.zeros((tp))
    for i in range(0,tp,1):
        acu=0
        for j in range(0,n,1):
            acu=acu+((PI[i,j]*PI[i,j])-A*np.cos(2*np.pi*PI[i,j]))
        fit[i]=acu+A*n    
    return fit

# Función de torneo
def torneo(PI,tp,n,fit):
    padres=np.zeros((tp,n))
    for i in range(0,tp,1):
        a=np.random.randint(tp)
        b=np.random.randint(tp)
        if fit[a]<fit[b]:
            padres[i,:]=PI[a,:]
        else:
            padres[i,:]=PI[b,:]
    return padres

# Función cruzamiento
def cruza(padres,tp,n):
    hijos=np.zeros((tp,n))
    for i in range(0,(tp),2):
        b=np.random.rand(1)
        if b <= pc:
            c=np.random.randint(1,n-1)
            # Hijo1 cabeza padre1 cola padre2
            hijos[i,0:c+1]=padres[i,0:c+1]
            hijos[i,c:n+1]=padres[i+1,c:n+1]
            # Hijo2 cabeza padre2
            hijos[i+1,0:c+1]=padres[i+1,0:c+1]
            hijos[i+1,c:n+1]=padres[i,c:n+1]
        else:
            hijos[i,:]=padres[i,:]
            hijos[i+1,:]=padres[i+1,:]
    return hijos

# Función mutación
def mut(hijos,tp,n,pm,a,b):
    hijosm=hijos.copy()
    for i in range(0,tp):
        b=np.random.rand(1)
        if b<= pm:
            c=np.random.randint(n)
            d=a+(b-a)*np.random.rand()
            hijosm[i,c]=d
    return hijosm

# Función actualizar óptimo
def optimo(PI,fit,opt,optfit,tp,n):
    for i in range(0,tp):
        if fit[i]<optfit:
            optfit=fit[i]
            opt[0,:]=PI[i,:]
    return opt,optfit

# Programa Principal
n = int(input('Introduce n: '))
tp = int(input('Introduce el tamaño de la población: '))
iter = int(input('Introduce el número de iteraciones: '))
pc = float(input('Introduce la probabilidad de cruzamiento: '))
pm = float(input('Introduce la probabilidad de mutación: '))
A=10
a=-5.12
b=5.12
opt=np.zeros((1,n))
optfit=10000000000000000000000
PI=PobIni(tp,n,a,b)
# Iteraciones
for i in range(0,iter+1):
    # Evaluación o fitness
    fit=fitness(PI,tp,n,A)
    # Actualización del óptimo
    (opt,optfit)=optimo(PI,fit,opt,optfit,tp,n)
    # Selección por torneo
    padres=torneo(PI,tp,n,fit)
    # Cruzamiento
    hijos=cruza(padres,tp,n)
    # Evaluación o fitness
    fit=fitness(hijos,tp,n,A)
    # Actualización del óptimo
    (opt,optfit)=optimo(hijos,fit,opt,optfit,tp,n)
    # Mutación
    hijosm=mut(hijos,tp,n,pm,a,b)
    # Reemplazo de la población
    PI=hijosm.copy()
    print("Iteración: ", i)
print('opt')
print(opt)
print('optfit')
print(optfit)