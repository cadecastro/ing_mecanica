#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 20 11:44:28 2021
COLUMNAS Y PANDEO
@author: cadecastro.com
"""
import numpy as np
import matplotlib.pyplot as plt
P=float(input('Carga aplicada P[N] ='))
L=float(input('Longitud columna L[mm] ='))
I=float(input('Momento de inercia menor I[mm⁴] ='))
A=float(input('Área transversal A[mm²] ='))
E=float(input('Módulo elasticidad E[MPa] ='))
Sy=float(input('Esfuerzo de fluencia Sy[MPa] ='))
print('Extremos: 1-articulados, 2-empotrados, 3-articulado-empotrado 4-empotrado-libre')
tipo=str(input('Tipo extremos ='))
#SOLUCIÓN:
#Longitud efectiva pandeo:
if tipo=='1' or tipo=='articulados':
    L=L
elif tipo=='2' or tipo=='empotrados':
    L=0.5*L
elif tipo=='3' or tipo=='articulado-empotrado':
    L=0.707*L
elif tipo=='4' or tipo=='empotrado-libre':
    L=2*L
else:
    print('Opción no válida, se tomará como extremos articulados.')
#SOLUCIÓN:
#Parámetros geométricos:
alpha=Sy/(np.pi*np.pi*E) #Constante Rankine
r=np.sqrt(I/A) #[mm] Radio de giro
S_r=L/r #Razón de esbeltez
#Esfuerzos:
sigma_apl=P/A #[Mpa] Esf. aplicado
sigma_euler=np.pi*np.pi*E/(S_r*S_r) #[Mpa] Esf. crítico Euler
sigma_rankine=Sy/(1+alpha*S_r*S_r) #[Mpa] Esf. crítico Rankine
print('RESULTADOS:')
print('Radio de giro =',np.format_float_positional(r,precision=2),'mm')
print('Razón de esbeltez =',np.format_float_positional(S_r,precision=2))
print('Esfuerzo crítico Euler =',np.format_float_positional(sigma_euler,precision=2),'MPa')
print('Esfuerzo crítico Rankine =',np.format_float_positional(sigma_rankine,precision=2),'MPa')
print('Factor seguridad Euler =',np.format_float_positional(sigma_euler/sigma_apl,precision=2))
print('Factor seguridad Rankine =',np.format_float_positional(sigma_rankine/sigma_apl,precision=2))
#Gráficas:
La=np.linspace(L/10,5*L,200)
Sr=np.linspace(S_r/10,5*S_r,200)
sigmaeuler=np.pi*np.pi*E/(Sr*Sr) #[Mpa] Esf. crítico Euler
for i in range(0,200):
    if sigmaeuler[i]>Sy:
        sigmaeuler[i]=Sy
sigmarankine=Sy/(1+alpha*Sr*Sr) #[Mpa] Esf. crítico Rankine
plt.figure(1)
plt.plot(Sr,sigmaeuler,'r')
plt.plot(Sr,sigmarankine,'b')
plt.plot(S_r,sigma_apl,'ko')
plt.title('Análisis pandeo',loc='left')
plt.title('cadecastro.com',loc='right')
plt.xlabel('Razón de esbeltez')
plt.ylabel('Esfuerzo [MPa]')
plt.grid(True,'both','both')
plt.legend(['Euler','Rankine','Aplicado'])

plt.figure(2)
plt.plot(La,sigmaeuler*A,'r')
plt.plot(La,sigmarankine*A,'b')
plt.plot(L,P,'ko')
plt.title('Análisis pandeo',loc='left')
plt.title('cadecastro.com',loc='right')
plt.xlabel('Longitud efectiva [mm]')
plt.ylabel('Carga axial [N]')
plt.grid(True,'both','both')
plt.legend(['Euler','Rankine','Aplicado'])