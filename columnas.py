def Pandeo_Columnas_SI( P , L, I , A , E , Sy, tipo ):

  print('Extremos: 1-articulados, 2-empotrados, 3-articulado-empotrado 4-empotrado-libre')
  #tipo=str(input('Tipo extremos ='))
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

  print(' ')
  print('RESULTADOS:')
  print('-------------------------------------------------------')
  print('Propiedades de la columna:')
  print('')
  print('           Radio de giro =',np.format_float_positional(r,precision=2),'mm')
  print('       Razón de esbeltez =',np.format_float_positional(S_r,precision=2))
  print('')
  print('-------------------------------------------------------')
  print('Esfuerzos:')
  print('')
  print('       Esfuerzo aplicado =',np.format_float_positional(sigma_apl,precision=2),'MPa')
  print('  Esfuerzo crítico Euler =',np.format_float_positional(sigma_euler,precision=2),'MPa')
  print('Esfuerzo crítico Rankine =',np.format_float_positional(sigma_rankine,precision=2),'MPa')
  print('')
  print('-------------------------------------------------------')
  print('Factores seguridad:')
  print('')
  print('  Factor seguridad Euler =',np.format_float_positional(sigma_euler/sigma_apl,precision=2))
  print('Factor seguridad Rankine =',np.format_float_positional(sigma_rankine/sigma_apl,precision=2))
  print('-------------------------------------------------------')
  print(' ')

  #Gráficas:
  La=np.linspace(L/10,1.5*L,200)
  Sr=np.linspace(S_r/10,1.5*S_r,200)
  sigmaeuler=np.pi*np.pi*E/(Sr*Sr) #[Mpa] Esf. crítico Euler
  for i in range(0,200):
      if sigmaeuler[i]>Sy:
          sigmaeuler[i]=Sy
  sigmarankine=Sy/(1+alpha*Sr*Sr) #[Mpa] Esf. crítico Rankine
  plt.figure(1 , figsize=(8,10) )
  plt.subplot(211)
  plt.fill_between(x=Sr ,y1=sigmaeuler , y2=sigmarankine , color='navy')
  plt.fill_between(x=Sr, y1=sigmarankine , y2=0 ,color='lime')
  plt.scatter( [S_r] , [sigma_apl] , marker='o', color='red')
  plt.title('Análisis esfuerzos de pandeo',loc='left', weight='bold', size=12)
  plt.xlabel('Razón de esbeltez')
  plt.ylabel('Esfuerzo [MPa]')
  plt.grid(True,'both','both')
  plt.legend(['Zona segura Euler','Zona segura Rankine','Aplicado'] , framealpha=0, prop={'size':8} , loc=1 )
  plt.xlim( 0 , None )
  plt.ylim( 0, None )

  plt.subplot(212)
  plt.fill_between( x=La , y1=sigmarankine*A, y2=0 , color='lime')
  plt.fill_between( x=La , y1=sigmaeuler*A , y2=sigmarankine*A ,color='navy')
  plt.scatter( [L] , [P] , marker='o', color='red' )
  plt.title('Análisis fuerzas de pandeo',loc='left', weight='bold', size=12)
  plt.title('cadecastro.com',style='italic', loc='right', color='navy', size=8)
  plt.xlabel('Longitud efectiva [mm]')
  plt.ylabel('Carga axial [N]')
  plt.grid(True,'both','both')
  plt.xlim( 0 , None )
  plt.ylim( 0, None )

print('Función análisis columnas definida.')
