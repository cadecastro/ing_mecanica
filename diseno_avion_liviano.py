#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 13 20:01:00 2021
DISEÑO Y DESEMPEÑO AVIÓN LIVIANO
@author: cadecastro.com
"""
import numpy as np
import matplotlib.pyplot as plt
#MASA DEL AVIÓN:
m=250#[kg] Masa del avión + Ocupantes + Carga

#DIMENSIONES ALA:
b=9.00#[m] Envergadura
c=1.74#[m] Cuerda

#MOTOR Y HÉLICE:
P_motor=46.9#[HP] Potencia nominal motor
eta_p=0.60# Eficiencia hélice
Tipo=2# Opciones: 1=combustión , 2=eléctrico

#PERFIL ALAR:
CL_max_perfil=2.353# Coef. Sustentación máximo
alpha1=15.0#[grados] Ángulo @ CL_max
L_D_max_perfil=17.8# (L/D)_max del perfil
CL1_perfil=1.479# Coef. Sustentación @ L/D_max perfil
alpha0=-9.0#[grados] Ángulo @ CL=0

#ARRASTRE PARÁSITO:
CD0_f=0.032# Coef. arrastre por fricción
CD0_otros=0.008# Coef. arrastre por otros factores

#CONDICIONES AMBIENTALES:
z=0#[m] Elevación sobre el nivel del mar
g=9.81#[m/s^2]
Tmar=30#[ºC] Temperatura al nivel del mar
Patm0=100000#[Pa] Presión atmosférica al nivel del mar
gradT=0.006562#[ºC/m] Gradiente de temperatura

#CONDICIONES SOLUCIONES NUMÉRICAS:
itmax=int(100)# Iteraciones máximas
e_max=1e-4# Residual relativo máximo
N=int(50)# Puntos a graficar

#--------------------------------------------------------------------------
#SOLUCIÓN:
#--------------------------------------------------------------------------
#Parámetros aerodinámicos perfil:
CD0_perfil=CL1_perfil/(2*L_D_max_perfil)#
K_perfil=1/(2*CL1_perfil*L_D_max_perfil)#
a0=CL_max_perfil/(alpha1-alpha0)*180/np.pi#[CL/rad] Pendiente curva CL perfil

#Parámetros aerodinámicos avión:
S=b*c#[m^2] Área alar
W=m*g#[N] Peso del avión
W_S=W/S#[N/m^2] Carga alar
AR=b/c# Razón de aspecto del ala
epsilon=1.78*(1-0.045*np.power(AR,0.68))-0.64# Factor de arrastre inducido
if epsilon>0.75:
    epsilon=0.75#
CL_max=0.85*CL_max_perfil# Coef. sustentación máximo avión
CD0=CD0_f+CD0_otros# Coef. arrastre parásito avión
K=1/(np.pi*epsilon*AR)+K_perfil#
L_D_max=0.5/np.sqrt(CD0*K)# (L/D)_max del avión
a=a0/(a0/(np.pi*AR)+np.sqrt(1+np.power(a0/(np.pi*AR),2)))#[CL/rad] Piente curva CL perfil
a0=a0*np.pi/180#[CL/grados] Piente curva CL perfil
a=a*np.pi/180#[CL/grados] Piente curva CL avión

#Curvas aerodinámicas:
Na=int(100)
alpha=np.linspace(alpha0,alpha1,Na)#
CL_perfil=np.zeros(Na)#
CD_perfil=np.zeros(Na)#
CL=np.zeros(Na)#
CD=np.zeros(Na)#
L_D_perfil=np.zeros(Na)#
L_D=np.zeros(Na)#
for i in range(0,Na):
    CL_perfil[i]=a0*(alpha[i]-alpha0)#
    CD_perfil[i]=CD0_perfil+K_perfil*CL_perfil[i]*CL_perfil[i]#
    L_D_perfil[i]=CL_perfil[i]/CD_perfil[i]#
    CL[i]=a*(alpha[i]-alpha0)#
    CD[i]=CD0+K*CL[i]*CL[i]#
    L_D[i]=CL[i]/CD[i]#

#Parámetros ambientales:
Tamb=Tmar-gradT*z+273.15#[K] Temperatura a la elevación dada
Patm=Patm0*np.power(1-gradT*z/(Tmar+273.15),g/(287*gradT))#[Pa] Presión atmosférica local
rho=Patm/(287*Tamb)#[kg/m^3] Densidad del aire local

#Potencia printonible motor:
if Tipo==2:
    Pdisp=eta_p*P_motor*746#[W] Potencia printonible motor eléctrico
else:
    Pdisp=eta_p*P_motor*746*(1.132*rho/1.225-0.132)#[W] Potencia printonible motor combustión


#Velocidades importantes:
V_stall=np.sqrt(2*W_S/(rho*CL_max))#[m/s] Velocidad de pérdida
Vy=np.sqrt(2*W_S/rho*np.sqrt(K/(3*CD0)))#[m/s] Velocidad máxima tasa de ascenso
V_GRmax=np.sqrt(2*W_S/rho*np.sqrt(K/CD0))#[m/s] Velocidad mejor planeo

#Cálculo velocidad máxima:
Vmax1=V_stall*1.8#Suposición inicial
cdc=1#
it=0#
while cdc==1:
   Vmax=Vmax1-(W*Vmax1*(0.5*rho*Vmax1*Vmax1*CD0/W_S+2*W_S*K/(rho*Vmax1*Vmax1))-Pdisp)/(W*(1.5*rho*Vmax1*Vmax1*CD0/W_S-2*W_S*K/(rho*Vmax1*Vmax1)))#[m/s] Vel. Máx.
   error=abs(Vmax-Vmax1)/Vmax1#
   it=it+1#
   if it>=itmax or error<=e_max:
       cdc=0#
   else:
       Vmax1=Vmax#

#Potencia y empuje mínimos requeridos:
Preq_min=1.755*W*np.sqrt(2*W_S/rho*np.sqrt(CD0*K*K*K))#[W] Potencia mínima requerida @ Vy
Treq_min=2*W*np.sqrt(CD0*K)#[N] Empuje mínimo requerido @ V_GRmax

#Máxima tasa de ascenso:
VS_max=(Pdisp-Preq_min)/W#[m/s] Máxima tasa de ascenso

#Cálculo curvas desempeño:
V=np.linspace(V_stall,Vmax,N)#
Treq=np.zeros(N)#
Preq=np.zeros(N)#
Pdispo=np.zeros(N)#
VS0=np.zeros(N)#
VS1=np.zeros(N)#
VS2=np.zeros(N)#
VS3=np.zeros(N)#
VS4=np.zeros(N)#
for i in range(0,N):
    Treq[i]=W*(rho*CD0*V[i]*V[i]/(2*W_S)+2*W_S*K/(rho*V[i]*V[i]))#[N]
    Preq[i]=V[i]*Treq[i]#[W]
    Pdispo[i]=Pdisp#[W]
    VS0[i]=(0-Preq[i])/W#[m/s]
    VS1[i]=(Preq_min-Preq[i])/W#[m/s]
    VS2[i]=(0.50*Pdisp-Preq[i])/W#[m/s]
    VS3[i]=(0.75*Pdisp-Preq[i])/W#[m/s]
    VS4[i]=(Pdisp-Preq[i])/W#[m/s]


#Viraje:
theta=np.linspace(0,60,61)
Nt=int(61)
n=np.zeros(Nt)#
Vstallv=np.zeros(Nt)#
for i in range(0,Nt):
   n[i]=1/np.cos(theta[i]*np.pi/180)# Factor de carga viraje
   Vstallv[i]=3600/1609*V_stall*np.sqrt(n[i])#[mph] Vel. pérdida en el viraje

#Conversión de unidades:
Preq=Preq*1/746#[HP]
Pdispo=Pdispo*1/746#[HP]
V=V*3600/1609#[mph]
VS0=VS0*60/0.3048#[ft/min]
VS1=VS1*60/0.3048#[ft/min]
VS2=VS2*60/0.3048#[ft/min]
VS3=VS3*60/0.3048#[ft/min]
VS4=VS4*60/0.3048#[ft/min]

#Salidas del programa:
print('Velocidad de pérdida=',np.format_float_positional(V_stall*3600/1609,precision=1),'mph')#
print('Velocidad de mejor tasa de ascenso Vy=',np.format_float_positional(Vy*3600/1609,precision=1),'mph')#
print('Velocidad de mejor tasa de planeo=',np.format_float_positional(V_GRmax*3600/1609,precision=1),'mph')#
print('Velocidad máxima=',np.format_float_positional(Vmax*3600/1609,precision=1),'mph')#
print('Mayor tasa de ascenso=',np.format_float_positional(VS_max*60/0.3048,precision=1),'ft/min')#
print('Mejor tasa de planeo =',np.format_float_positional(L_D_max,precision=1))#
print('Potencia mínima requerida=',np.format_float_positional(Preq_min*1/746,precision=1),'HP')#
print('Empuje mínimo requerido=',np.format_float_positional(Treq_min,precision=1),'N')#

#Gráficas:
plt.figure(1)
plt.subplot(141)
plt.plot(alpha,CL,'b')
plt.plot(alpha,CL_perfil,'r-.')
plt.grid(True,'both','both')
plt.title('Coef. sustentación vs. AOA')
plt.xlabel('Ángulo de ataque [grados]')
plt.ylabel('C_L')
plt.legend(['Avión','Perfil'])

plt.subplot(142)
plt.plot(alpha,CD,'b')
plt.plot(alpha,CD_perfil,'r-.')
plt.grid(True,'both','both')
plt.title('Coef. arrastre vs. AOA')
plt.xlabel('Ángulo de ataque [grados]')
plt.ylabel('C_D')
plt.legend(['Avión','Perfil'])

plt.subplot(143)
plt.plot(CL,L_D,'b')
plt.plot(CL_perfil,L_D_perfil,'r-.')
plt.grid(True,'both','both')
plt.title('L/D vs. C_L')
plt.xlabel('C_L')
plt.ylabel('L/D')
plt.legend(['Avión','Perfil'])

plt.subplot(144)
plt.plot(CD,CL,'b')
plt.plot(CD_perfil,CL_perfil,'r-.')
plt.grid(True,'both','both')
plt.title('Polar de arrastre')
plt.xlabel('C_D')
plt.ylabel('C_L')
plt.legend(['Avión','Perfil'])

plt.figure(2)
plt.subplot(211)
plt.plot(V,Preq,'b')
plt.plot(V,Pdispo,'r-.')
plt.grid(True,'both','both')
plt.title('Vuelo recto y nivelado')
plt.ylabel('Potencia [HP]')
plt.legend(['Requerida','Disponible'])

plt.subplot(212)
plt.plot(V,Treq,'b')
plt.grid(True,'both','both')
plt.xlabel('Velocidad [mph] - cadecastro.com')
plt.ylabel('Empuje requerido [N]')

plt.figure(3)
plt.plot(V,VS4,'r-')
plt.plot(V,VS3,'b--')
plt.plot(V,VS2,'r-.')
plt.plot(V,VS1,'b:')
plt.plot(V,VS0,'r-d')
plt.grid(True,'both','both')
plt.title('Curvas de desempeño')
plt.xlabel('Velocidad [mph]')
plt.ylabel('Velocidad vertical [pies/min]')
plt.legend(['100 % Potencia','75 % Potencia','50 % Potencia','Mínima requerida','Planeo'])

plt.figure(4)
plt.subplot(211)
plt.plot(theta,Vstallv,'b')
plt.grid(True,'both','both')
plt.title('V_stall y factor de carga en viraje - cadecastro.com')
plt.ylabel('Velocidad de pérdida [mph]')

plt.subplot(212)
plt.plot(theta,n,'b')
plt.grid(True,'both','both')
plt.xlabel('Ángulo de alabeo [grados]')
plt.ylabel('Factor de carga')