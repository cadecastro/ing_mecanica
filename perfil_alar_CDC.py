#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 13 11:46:59 2021

@author: cadecastro.com
"""
import numpy as np
import matplotlib.pyplot as plt
#%INGRESO DE PARÁMETROS:
t=0.15#;% Espesor como fracción de la cuerda al 25 % de la misma
c=1740#;%[mm] Cuerda del ala
N=int(1000)#Número de puntos gráfica
#%CURVAS DEL PERFIL:
R1=c*(16*t*t+1)/(32*t)
R2=c*(16*t*t+9)/(32*t)
x=np.linspace(0,c,N)
y_inf=np.zeros(N)
y_sup=np.zeros(N)
for i in range(0,N):
    if x[i]<=0.25*c:
        y_sup[i]=np.sqrt(R1*R1-(0.25*c-x[i])*(0.25*c-x[i]))-R1+t*c
    else:
        y_sup[i]=np.sqrt(R2*R2-(0.25*c-x[i])*(0.25*c-x[i]))-R2+t*c
#Salidas:
print('Cuerda [mm] = ',c)
print('Radio R1 [mm] = ',R1)
print('Radio R2 [mm] = ',R2)
#Gráfica:
plt.plot(x,y_inf,'k-')
plt.plot(x,y_sup,'k-')
plt.grid(True,'both','both')
plt.axis('equal')
plt.title('Gráfica del perfil alar')
plt.xlabel('x [mm]')
plt.ylabel('y [mm]')