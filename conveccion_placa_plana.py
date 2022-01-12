#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CONVECCIÓN PLACA PLANA EN AIRE
Created on Thu Nov  4 21:59:46 2021
@author: cadecastro.com
"""
import numpy as np
L=float(input('Longitud placa L [m]=')) #[m] Longitud de la placa
b=float(input('Ancho placa b [m]=')) #Ancho de la placa [m]
T_inf=float(input('Temp. del aire T_inf [C]=')) #[C] Temperatura alrededores
Ts=float(input('Temp. superficial de la placa Ts [C]=')) #[C] Temperatura superficial
V_inf=float(input('Velocidad del aire V_inf [m/s]='))#[m/s] Velocidad aire
Pref=float(input('Presión de referencia P_ref [kPa]=')) #[kPa] Presión de referencia

n=3.5 #Factor Nusselt convección mixta
g=9.81 #[m/s²] Gravedad
#SOLUCIÓN:
A=L*b #[m²] Área placa
P=2*(L+b) #[m] Perímetro placa
Lc=A/P #[m] Long. característica convección natural
Tf=0.5*(T_inf+Ts)+273.15 #[K] Temperatura de película
beta=1/Tf #[1/K] Coef. expansión
rho=Pref/(0.287*Tf) #[kg/m³] Densidad del aire
Cp=1003.62 #[J/kg*K] Calor específico aire
mu=1.716e-5*np.power(Tf/273.15,1.5)*384.15/(Tf+111) #[Pa*s] Viscosidad dinámica
k=0.02414*np.power(Tf/273.15,1.5)*467.15/(Tf+194) #[W/m*K] Conductividad térmica
alpha=k/(rho*Cp) #[m²/s] Difusividad térmica
nu=mu/rho #[m²/s] Viscosidad cinemática
Pr=nu/alpha #Número de Prandtl
Re=V_inf*L/nu #Reynolds
Gr=g*beta*abs(Ts-T_inf)*np.power(Lc,3)/np.power(nu,2) #Grashof
Ra=Gr*Pr #Número de Rayleigh
if Re!=0:
    Ri=Gr/(Re*Re) #Número de Richardson
else:
    Ri=100
    print('No hay convección forzada')
if Ri>=10:
    print('Convección natural dominante')
    if Ra<=1e7:
        print('Convección natural laminar')
        Nu=0.54*np.power(Ra,0.25)
        h=k*Nu/Lc #[W/m²*K]
        if Ra<1e4:
            print('Ra<10^4 , correlación puede fallar')
        if Pr<0.7:
            print('Pr<0.7 , correlación puede fallar')
    else:
        print('Convección natural turbulenta')
        Nu=0.15*np.power(Ra,1/3)
        h=k*Nu/Lc #[W/m²*K]
        if Ra>1e11:
            print('Ra>10^11 , correlación puede fallar')
elif Ri<=0.1:
    print('Convección forzada dominante')
    if Re<5e5:
        print('Convección forzada laminar')
        Nu=0.664*np.power(Re,0.5)*np.power(Pr,1/3)
        h=k*Nu/L #[W/m²*K]
        if Pr<0.6 or Pr>50:
            print('Pr fuera de límites, correlación puede fallar')
    else:
        print('Convección forzada turbulenta')
        Nu=(0.037*np.power(Re,4/5)-871)*np.power(Pr,1/3)
        h=k*Nu/L #[W/m²*K]
else:
    print('Convección natural y forzada mixtas')
    if Ra<=1e7:
        print('Convección natural laminar')
        Nu_nat=0.54*np.power(Ra,0.25)
        h_nat=k*Nu_nat/Lc #[W/m²*K]
        if Ra<1e4:
            print('Ra<10^4 , correlación puede fallar')
        if Pr<0.7:
            print('Pr<0.7 , correlación puede fallar')
    else:
        print('Convección natural turbulenta')
        Nu_nat=0.15*np.power(Ra,1/3)
        if Ra>1e11:
            print('Ra>10^11 , correlación puede fallar')
    if Re<5e5:
        print('Convección forzada laminar')
        Nu_f=0.664*np.power(Re,0.5)*np.power(Pr,1/3)
        if Pr<0.6 or Pr>50:
            print('Pr fuera de límites, correlación puede fallar')
    else:
        print('Convección forzada turbulenta')
        Nu_f=(0.037*np.power(Re,4/5)-871)*np.power(Pr,1/3)
    Nu=np.power(np.power(Nu_nat,n)+np.power(Nu_f,n),1/n)
    h=k*Nu/Lc
q=h*(Ts-T_inf) #[W/m²]
Q=q*A #[W]
#Resultados:
print('RESULTADOS:')
print('Pr=',np.format_float_positional(Pr,precision=3))
print('Re=',np.format_float_scientific(Re,precision=3))
print('Gr=',np.format_float_scientific(Gr,precision=3))
print('Ra=',np.format_float_scientific(Ra,precision=3))
if Re==0:
    print('No aplica Ri')
else:
    print('Ri=',np.format_float_positional(Ri,precision=3))
print('Nu=',np.format_float_positional(Nu,precision=3))
print('Coef. convección h=',np.format_float_positional(h,precision=3),'W/m²*K')
print('Flujo de calor q"=',np.format_float_positional(q,precision=3),'W/m²')
print('Tasa de calor dQ/dt=',np.format_float_positional(Q,precision=3),'W')