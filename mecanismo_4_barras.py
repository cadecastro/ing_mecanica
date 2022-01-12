#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ANÁLISIS MECANISMO DE CUATRO BARRAS
Created on Sun Nov 14 09:14:08 2021
@author: cadecastro.com
"""
import numpy as np
import matplotlib.pyplot as plt
#Medidas de barras del mecanismo:
L2=0.300 #[m] Long. barra 2 (entrada)
L3=0.500 #[m] Long. barra 3
L4=0.150 #[m] Long. barra 4
Z1=0.000 #[m] Long. de A a P
beta=0.00 #[grados] Ángulo de AP respecto a barra 3
#Coordenadas posición apoyo O2 (O1 se toma en el origen)
O2x=0.583 #[m]
O2y=0.050 #[m]
#Ángulo de entrada
theta2i=45 #[grados] Ángulo entrada inicial
theta2f=135 #[grados] Ángulo entrada final
N=int(91) #Número de posiciones a analizar
w2=6 #[rad/s] Velocidad angular barra 2
#Parámetros solución numérica posición:
itmax=100 #Iteraciones máximas por ángulo
tolerancia=1e-6 #[rad] Tolerancia error absoluto
f=0.5 #Factor suposición inicial

#SOLUCIÓN:
theta2=np.linspace(theta2i,theta2f,N)*np.pi/180 #[rad]
theta3=np.zeros(N)
theta4=np.zeros(N)
beta=beta*np.pi/180 #[rad]
for i in range(0,N):
    it=0
    residual=1
    cdc=1
    theta31=f*theta2[i] #Suposición inicial theta3
    while cdc==1:
        it=it+1
        num1=L2*np.sin(theta2[i])+L3*np.sin(theta31)-O2y
        den1=L2*np.cos(theta2[i])+L3*np.cos(theta31)-O2x
        if den1==0:
            theta41=np.pi/2
        else:
            if den1<0 and num1>0:
                theta41=np.pi+np.arctan(num1/den1)
            elif den1<0 and num1<0:
                theta41=np.pi+np.arctan(num1/den1)
            elif den1>0:
                theta41=np.arctan(num1/den1)
        num1=O2y-L2*np.sin(theta2[i])+L4*np.sin(theta41)
        den1=O2x-L2*np.cos(theta2[i])+L4*np.cos(theta41)
        if den1==0:
            theta32=np.pi/2
        else:
            if den1<0 and num1>0:
                theta32=np.pi+np.arctan(num1/den1)
            elif den1<0 and num1<0:
                theta32=np.pi+np.arctan(num1/den1)
            elif den1>0:
                theta32=np.arctan(num1/den1)
        residual=abs(theta32-theta31)
        if it>=itmax or residual<=tolerancia:
            cdc=0
            theta3[i]=theta32
            theta4[i]=theta41
        else:
            theta31=theta32
    if it==itmax:
        print('Alerta: máximas iteraciones alcanzadas para theta2=',np.format_float_positional(theta2[i]*180/np.pi,precision=0),'grados')
        print('Residual absoluto=',np.format_float_scientific(residual,precision=3),'rad')
#Coordenadas punto A:
Ax=L2*np.cos(theta2)
Ay=L2*np.sin(theta2)
#Coordenadas punto B:
Bx=O2x+L4*np.cos(theta4)
By=O2y+L4*np.sin(theta4)
#Coordenadas punto de salida P:
Px=L2*np.cos(theta2)+Z1*np.cos(theta3+beta)
Py=L2*np.sin(theta2)+Z1*np.sin(theta3+beta)
#Análisis de velocidad:
vAx=-w2*L2*np.sin(theta2)
vAy=w2*L2*np.cos(theta2)
w3=(vAx*np.cos(theta4)+vAy*np.sin(theta4))/(L3*np.sin(theta3-theta4))
w4=(vAx*np.cos(theta3)+vAy*np.sin(theta3))/(L4*np.sin(theta3-theta4))
vBx=-w4*L4*np.sin(theta4)
vBy=w4*L4*np.cos(theta4)
vPx=vAx-w3*Z1*np.sin(theta3+beta)
vPy=vAy+w3*Z1*np.cos(theta3+beta)
#Análisis de aceleración:
aAx=-w2*vAy
aAy=w2*vAx
a1=aAx-w3*(vBy-vAy)+w4*vBy
a2=aAy+w3*(vBx-vAx)-w4*vBx
alpha3=(a1*np.cos(theta4)+a2*np.sin(theta4))/(L3*np.sin(theta3-theta4))
alpha4=(a1*np.cos(theta3)+a2*np.sin(theta3))/(L4*np.sin(theta3-theta4))
aBx=-alpha4*L4*np.sin(theta4)-w4*vBy
aBy=alpha4*L4*np.cos(theta4)+w4*vBx
aPx=aAx-alpha3*Z1*np.sin(theta3+beta)-w3*(vPy-vAy)
aPy=aAy+alpha3*Z1*np.cos(theta3+beta)+w3*(vPx-vAx)
#Gráficas de posición:
plt.figure(1)
plt.axis('equal')
plt.plot(Ax,Ay,'b--')
plt.plot(Px,Py,'r--')
plt.plot(Bx,By,'y--')
plt.plot(0,0,'ko')
plt.plot(O2x,O2y,'go')
plt.grid(True,'both','both')
plt.legend(['Tray. Punto A','Tray. Punto P','Tray. Punto B','Apoyo 1','Apoyo 2'])
#Mecanismo en la posición inicial:
plt.plot([0,Ax[0]],[0,Ay[0]],'k')
plt.plot([Ax[0],Px[0]],[Ay[0],Py[0]],'k')
plt.plot([Px[0],Bx[0]],[Py[0],By[0]],'k')
plt.plot([Ax[0],Bx[0]],[Ay[0],By[0]],'k')
plt.plot([Bx[0],O2x],[By[0],O2y],'k')
#Mecanismo en la posición media:
n=int(np.round(N/2,decimals=0))
plt.plot([0,Ax[n]],[0,Ay[n]],'c')
plt.plot([Ax[n],Px[n]],[Ay[n],Py[n]],'c')
plt.plot([Px[n],Bx[n]],[Py[n],By[n]],'c')
plt.plot([Ax[n],Bx[n]],[Ay[n],By[n]],'c')
plt.plot([Bx[n],O2x],[By[n],O2y],'c')
#Mecanismo en la posición final:
plt.plot([0,Ax[N-1]],[0,Ay[N-1]],'m')
plt.plot([Ax[N-1],Px[N-1]],[Ay[N-1],Py[N-1]],'m')
plt.plot([Px[N-1],Bx[N-1]],[Py[N-1],By[N-1]],'m')
plt.plot([Ax[N-1],Bx[N-1]],[Ay[N-1],By[N-1]],'m')
plt.plot([Bx[N-1],O2x],[By[N-1],O2y],'m')
plt.title('Posiciones mecanismo de 4 barras',loc='left')
plt.title('cadecastro.com',loc='right')
plt.xlabel('Negro: posición inicial - Cyan: pos. media - Magenta: pos. final')

plt.figure(2)
plt.plot(theta2*180/np.pi,theta3*180/np.pi,'b')
plt.plot(theta2*180/np.pi,theta4*180/np.pi,'r')
plt.grid(True,'both','both')
plt.title('Ángulos barras',loc='left')
plt.title('cadecastro.com',loc='right')
plt.legend(['theta3','theta4'])
plt.xlabel('Ángulo de entrada [grados]')
plt.ylabel('Ángulo de salida [grados]')
#Gráficas de velocidad:
plt.figure(3)
plt.plot(theta2*180/np.pi,vAx,'c-.')
plt.plot(theta2*180/np.pi,vAy,'m-.')
plt.plot(theta2*180/np.pi,vBx,'g--')
plt.plot(theta2*180/np.pi,vBy,'y--')
plt.plot(theta2*180/np.pi,vPx,'b')
plt.plot(theta2*180/np.pi,vPy,'r')
plt.grid(True,'both','both')
plt.title('Componentes velocidad puntos',loc='left')
plt.title('cadecastro.com',loc='right')
plt.legend(['vAx','vAy','vBx','vBy','vPx','vPy'])
plt.xlabel('Ángulo de entrada [grados]')
plt.ylabel('Velocidad [m/s]')

plt.figure(4)
plt.plot(theta2*180/np.pi,w3,'b')
plt.plot(theta2*180/np.pi,w4,'r')
plt.grid(True,'both','both')
plt.title('Velocidades angulares',loc='left')
plt.title('cadecastro.com',loc='right')
plt.legend(['w3','w4'])
plt.xlabel('Ángulo de entrada [grados]')
plt.ylabel('Velocidad angular [rad/s]')

#Gráficas de velocidad:
plt.figure(5)
plt.plot(theta2*180/np.pi,aAx,'c-.')
plt.plot(theta2*180/np.pi,aAy,'m-.')
plt.plot(theta2*180/np.pi,aBx,'g--')
plt.plot(theta2*180/np.pi,aBy,'y--')
plt.plot(theta2*180/np.pi,aPx,'b')
plt.plot(theta2*180/np.pi,aPy,'r')
plt.grid(True,'both','both')
plt.title('Componentes aceleración puntos',loc='left')
plt.title('cadecastro.com',loc='right')
plt.legend(['aAx','aAy','aBx','aBy','aPx','aPy'])
plt.xlabel('Ángulo de entrada [grados]')
plt.ylabel('Aceleración [m/s²]')

plt.figure(6)
plt.plot(theta2*180/np.pi,alpha3,'b')
plt.plot(theta2*180/np.pi,alpha4,'r')
plt.grid(True,'both','both')
plt.title('Aceleraciones angulares',loc='left')
plt.title('cadecastro.com',loc='right')
plt.legend(['alpha3','alpha4'])
plt.xlabel('Ángulo de entrada [grados]')
plt.ylabel('Aceleración angular [rad/s²]')