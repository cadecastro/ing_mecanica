#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 13 15:37:21 2021
@author: cadecastro
"""
#«Todos amamos a Python, y a su glorioso régimen»
import numpy as np
import matplotlib.pyplot as plt
v0=350#[m/s] Velocidad inicial
theta0=0#[grados] Ángulo de tiro
y0=1.50#[m] Altura inicial
tf=0.726#[s] Tiempo final
dt=0.0001#[s] Paso de tiempo
m=0.00065#[kg] Masa del cuerpo
S=np.pi/4*0.010*0.010#[m^2] Área superficial flujo
CD=0.5#Coeficiente de arrastre
CL=1.0#Coeficiente de sustentación
rho=1.225#[kg/m^3] Densidad del aire
g=9.81#[m/s^2] Magnitud aceleración gravedad
#SOLUCIÓN:
N=int(tf/dt+1)
theta0=theta0*np.pi/180#[rad]
t=np.linspace(0,tf,N)
x=np.zeros(N)#
y=np.zeros(N)#
vx=np.zeros(N)#
vy=np.zeros(N)#
v=np.zeros(N)#
L=np.zeros(N)#
D=np.zeros(N)#
theta=np.zeros(N)#
y[0]=y0#
vx[0]=v0*np.cos(theta0)#
vy[0]=v0*np.sin(theta0)#
v[0]=v0#
theta[0]=theta0#
L[0]=0.5*rho*v[0]*v[0]*S*CL#
D[0]=0.5*rho*v[0]*v[0]*S*CD#
for i in range(0,N-1):
    vx[i+1]=vx[i]-dt/m*(L[i]*np.sin(theta[i])+D[i]*np.cos(theta[i]))#
    vy[i+1]=vy[i]+dt/m*(L[i]*np.cos(theta[i])-D[i]*np.sin(theta[i]))-g*dt#
    x[i+1]=x[i]+dt*vx[i+1]#
    y[i+1]=y[i]+dt*vy[i+1]#
    if vx[i+1]!=0:
        theta[i+1]=np.arctan(vy[i+1]/vx[i+1])#
    else:
        if vy[i+1]>=0:
            theta[i+1]=np.pi/2#
        else:
            theta[i+1]=-np.pi/2#
    v[i+1]=np.sqrt(vx[i+1]*vx[i+1]+vy[i+1]*vy[i+1])#
    L[i+1]=0.5*rho*v[i+1]*v[i+1]*S*CL#
    D[i+1]=0.5*rho*v[i+1]*v[i+1]*S*CD#
theta=theta*180/np.pi#
if theta0==0:
    print('Disparo horizontal')
else:
    print('Disparo a',theta0*180/np.pi,'grados')
    print('Altura máxima =',np.max(y),'m')
print('Alcance horizontal =',x[N-1],'m')
print('Altura final =',y[N-1],'m')
print('Velocidad final =',v[N-1],'m/s')
print('Energía inicial =',0.5*m*v[0]*v[0],'J')
print('Energía final =',0.5*m*v[N-1]*v[N-1],'J')
print('Pérdida de energía =',(1-v[N-1]*v[N-1]/(v[0]*v[0]))*100,'%')
#Gráficas:
plt.figure(1)
plt.grid(True,'both','both')
#plt.axis('equal')
plt.title('Trayectoria')
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.plot(x,y)

#plt.figure(2)
#plt.plot(t,x,'b')
#plt.plot(t,y,'r')
#plt.grid(True,'both','both')
#plt.title('Posición')
#plt.xlabel('t [s]')
#plt.ylabel('Posición [m]')
#plt.legend(['x','y'])

#plt.figure(3)
#plt.plot(t,vx,'b')
#plt.plot(t,vy,'r')
#plt.plot(t,v,'g')
#plt.grid(True,'both','both')
#plt.title('Velocidad')
#plt.xlabel('t [s]')
#plt.ylabel('Velocidad [m/s]')
#plt.legend(['v_x','v_y','|v|'])

#plt.figure(4)
#plt.plot(t,theta,'b')
#plt.grid(True,'both','both')
#plt.title('Dirección')
#plt.xlabel('t [s]')
#plt.ylabel('\theta [grados]')
#plt.figure(5)
#plt.plot(t,L,'b')
#plt.plot(t,D,'r')
#plt.grid(True,'both','both')
#plt.title('Fuerzas del aire')
#plt.xlabel('t [s]')
#plt.ylabel('Fuerza [N]')
#plt.legend(['Sustentación','Arrastre'])

plt.figure(5)
plt.plot(x,v,'r')
plt.grid(True,'both','both')
plt.title('Velocidad vs. distancia horizontal')
plt.xlabel('x [m]')
plt.ylabel('Velocidad [m/s]')