#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 17 17:04:05 2021
@author: cadecastro.com
"""
import numpy as np
import matplotlib.pyplot as plt
#Estado de esfuerzos:
sigma_x=float(input('sigma_x [MPa] ='))
sigma_y=float(input('sigma_y [MPa] =')) #[MPa]
tau_xy=float(input('tau_xy [MPa] =')) #[MPa]
#Propiedades material:
print('¿Analizar falla? Dúctil=1, Frágil=2, No=otro valor')
tipo=str(input('Opción: '))
if tipo=='1':
    Sy=float(input('Sy [MPa] ='))
elif tipo=='2':
    Sut=float(input('Sut [MPa] ='))
    Suc=float(input('Suc [MPa] ='))
#Esfuerzos principales:
sigma_m=0.5*(sigma_x+sigma_y)
tau_max=np.sqrt(np.power(0.5*(sigma_x-sigma_y),2)+tau_xy*tau_xy)
sigma1=sigma_m+tau_max
sigma2=sigma_m-tau_max
theta_p=0.5*np.arctan(2*tau_xy/(sigma_x-sigma_y))
theta_s=0.5*np.arctan((-sigma_x+sigma_y)/(2*tau_xy))
print('Esfuerzos principales:')
print('sigma1=',np.format_float_positional(sigma1,precision=2),'MPa')
print('sigma2=',np.format_float_positional(sigma2,precision=2),'MPa')
print('tau_max=',np.format_float_positional(tau_max,precision=2),'MPa')
print('Orientación principales:')
print('theta_p=',np.format_float_positional(theta_p*180/np.pi,precision=2),'grados')
print('theta_s=',np.format_float_positional(theta_s*180/np.pi,precision=2),'grados')
#Círculo de Mohr:
sigma=np.linspace(sigma2,sigma1,100)
tau_sup=np.sqrt(tau_max*tau_max-np.power(sigma-sigma_m,2))
tau_inf=-np.sqrt(tau_max*tau_max-np.power(sigma-sigma_m,2))
plt.figure(1)
plt.plot(sigma1,0,'ro')
plt.plot(sigma2,0,'yo')
plt.plot(sigma,tau_sup,'b')
plt.plot(sigma,tau_inf,'b')
plt.plot(sigma_x,tau_xy,'ko')
plt.plot(sigma_y,-tau_xy,'ko')
plt.plot(sigma_m,0,'ko')
plt.plot([sigma_x,sigma_y],[tau_xy,-tau_xy],'k--')
plt.axis('equal')
plt.title('Círculo de Mohr',loc='left')
plt.title('cadecastro.com',loc='right')
plt.xlabel('sigma [MPa]')
plt.ylabel('tau [MPa]')
plt.grid(True,'both','both')
plt.legend(['sigma1','sigma2'])
#Análisis de falla:
if tipo=='1':
    sigma_vm=np.sqrt(sigma_x*sigma_x-sigma_x*sigma_y+sigma_y*sigma_y+3*tau_xy*tau_xy) #[MPa] Esfuerzo de von Mises
    tau_y=0.5*Sy
    n_vm=Sy/sigma_vm
    n_tresca=tau_y/tau_max
    sigma1vm=np.linspace(-Sy*np.sqrt(4/3),Sy*np.sqrt(4/3),100)
    sigma2vms=0.5*(sigma1vm+np.sqrt(sigma1vm*sigma1vm-4*(sigma1vm*sigma1vm-Sy*Sy)))
    sigma2vmi=0.5*(sigma1vm-np.sqrt(sigma1vm*sigma1vm-4*(sigma1vm*sigma1vm-Sy*Sy)))
    print('Análisis de falla dúctil:')
    print('Teoría de von Mises:')
    print('sigma_vm=',np.format_float_positional(sigma_vm,precision=2),'MPa')
    print('Factor de seguridad Von Mises=',np.format_float_positional(n_vm,precision=3))
    print('Teoría de Tresca (máximo cortante):')
    print('Tresca tau_y=',np.format_float_positional(tau_y,precision=2),'MPa')
    print('Factor de seguridad Tresca=',np.format_float_positional(n_tresca,precision=3))
    #Gráficas de locaciones de falla:
    plt.figure(2)
    plt.plot(sigma1,sigma2,'ro')
    plt.plot([-Sy,0,Sy,Sy,0,-Sy,-Sy],[0,Sy,Sy,0,-Sy,-Sy,0],'b')
    plt.plot(sigma1vm,sigma2vms,'g')
    plt.plot(sigma1vm,sigma2vmi,'g')
    plt.axis('equal')
    plt.title('Análisis de falla dúctil',loc='left')
    plt.title('cadecastro.com',loc='right')
    plt.xlabel('sigma1 [MPa]')
    plt.ylabel('sigma2 [MPa]')
    plt.grid(True,'both','both')
    plt.legend(['Estado de esfuerzos','Tresca','von Mises'])
elif tipo=='2':
    if sigma1>0 and sigma2<0:
        n=np.power(sigma1/Sut-sigma2/Suc,-1)
    elif sigma1<0:
        n=-Suc/sigma2
    elif sigma2>0:
        n=Sut/sigma1
    print('Análisis de falla frágil:')
    print('Teoría de Coulomb-Mohr:')
    print('Factor de seguridad=',np.format_float_positional(n,precision=3))
    #Gráficas de locaciones de falla:
    plt.figure(2)
    plt.plot(sigma1,sigma2,'ro')
    plt.plot([-Suc,0,Sut,Sut,0,-Suc,-Suc],[0,Sut,Sut,0,-Suc,-Suc,0],'b')
    plt.axis('equal')
    plt.title('Análisis de falla frágil',loc='left')
    plt.title('cadecastro.com',loc='right')
    plt.xlabel('sigma1 [MPa]')
    plt.ylabel('sigma2 [MPa]')
    plt.grid(True,'both','both')
    plt.legend(['Estado de esfuerzos','Coulomb-Mohr'])
else:
    print('No se analizó falla del material.')
