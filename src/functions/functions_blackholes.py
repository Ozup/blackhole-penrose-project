import numpy as np
import matplotlib.pyplot as plt
from . import functions_aux as fa

import importlib
importlib.reload(fa)

# Funciones de la sección 11.3 | Figura 11.1

def plotlightcone(r, C, col="lightcoral", mu = 1, linew = 0.6):
  """
  Grafica las lineas de mundo de fotones entrantes y saleintes en coordenadas de Schwarzschild.
  r: Distancias radiales que se desea graficar. (array)
  C: Constante arbitraría que desplaza hacia arriba o hacia abajo las curvas. (float)
  col: Color de la gráfica (Opcinal, está por defecto en lightcoral). (str)
  mu: Es la abreviación de GM (Opcional, trabajamos normalmente con mu=1). (float)
  """
  plt.plot(r, fa.ct(r, C, True, mu=mu), color = col, linewidth= linew)
  plt.plot(r, fa.ct(r, C, False, mu=mu), color = col, linestyle="--", linewidth= linew)
  return

# Funciones de la sección 11.4 | Figura 11.2

def time_11_4(r0, r, mu): # | Ecuación 11.4 | Capítulo:11 | Sección:11.4 | Página 253 |
  """
  Ecuación de linea de mundo de una partícula con masa entrante en un agujero negro en coordenadas de Schwarzschild.
  Retorna un array de tiempo correspondiente a los r elegidos.
  r0: Distancia radial donde t=0 y tao=0. (float)
  r: Distancias radiales que se desea graficar. (array)
  mu: Es la abreviación de GM (Opcional, trabajamos normalmente con mu=1). (float)
  """
  c = 1
  term_1 = 2/3 * (np.sqrt((r0**3)/(2*mu*c**2)) - np.sqrt((r**3)/(2*mu*c**2)))
  term_2 = 4*mu/c * (np.sqrt((r0)/(2*mu)) - np.sqrt((r)/(2*mu)))
  term_3 = 2*mu/c * np.log(abs( ((np.sqrt(r/(2*mu)) + 1)/(np.sqrt(r/(2*mu)) - 1)) * ((np.sqrt(r0/(2*mu)) - 1)/(np.sqrt(r0/(2*mu)) + 1)) ))
  return (term_1 + term_2 + term_3)*c/mu

def plot_tao(r0, tao, mu, off_setx=0.1, off_sety=0):
  """
  Grafica un punto con un texto al lado del tao correspondiente.
  r0: Distancia radial donde t=0 y tao=0. (float)
  tao: Tiempo propio del que se desea hallar el r y t correspondiente. (float)
  mu: Es la abreviación de GM (Opcional, trabajamos normalmente con mu=1). (float)
  off_setx: Valor que desea que se mueva el texto en x. (Opcional, por defecto está en 0.1) (float)
  off_sety: Valor que desea que se mueva el texto en y. (Opcional, por defecto está en 0) (float)
  """
  plt.plot(fa.time_tao(r0, tao, mu)[0], fa.time_tao(r0, tao, mu)[1], "o", markersize="3", color = "black")
  plt.text(fa.time_tao(r0, tao, mu)[0] + off_setx, fa.time_tao(r0, tao, mu)[1] + off_sety, r"$\tau= {}\mu/c$".format(tao) )
  return

def plotlightcone1(r, C1, C2, tao, r0, col="steelblue", mu = 1):
  """
  Grafica un cono al rededor de un punto tao de interés.
  r: Distancias radiales que se desea graficar. (array)
  C1: Constante de apoyo para graficar el cono. (float)
  C2: Constante de apoyo para graficar el cono. (float)
  tao: Tiempo propio del que se desea hallar el r y t correspondiente. (float)
  r0: Distancia radial donde t=0 y tao=0. (float)
  mu: Es la abreviación de GM (Opcional, trabajamos normalmente con mu=1). (float)
  """
  x, y, mu = fa.time_tao(r0, tao, mu)
  n = 0.5
  r = np.linspace(x-n, x+n, 100)
  plt.plot(r, fa.ct(r, C1, True), color = col, linestyle="--", linewidth= 0.7)
  plt.plot(r, fa.ct(r, C2, False), color = col, linestyle="--",  linewidth= 0.7)
  return

def Buscar_C(r0, tao, mu):
  """
  Busca el C adecuado para graficar el cono en el punto de tao de interés.
  Retorna los C1 y C2 auxiliares para la busqueda el cono en el tao de interés.
  r0: Distancia radial donde t=0 y tao=0. (float)
  tao: Tiempo propio del que se desea hallar el r y t correspondiente. (float)
  mu: Es la abreviación de GM (Opcional, trabajamos normalmente con mu=1). (float)
  """
  c = 1
  r, t, mu = fa.time_tao(r0, tao, mu)
  C1 = -r - 2*mu*np.log(r/(2*mu)-1) + c*t
  C2 = r + 2*mu*np.log(r/(2*mu)-1) + c*t
  return C1, C2

# Funciones de la sección 11.5 | Figura 11.3

def plotlightcone2(r, C, col="skyblue", mu = 1, zone2=True, linew=0.7):
  """
  Grafica las lineas de mundo de fotones entrantes y saleintes en coordenadas de Eddington-Filkelstein Avanzadas.
  r: Distancias radiales que se desea graficar. (array)
  C: Constante arbitraría que desplaza hacia arriba o hacia abajo las curvas. (float)
  col: Color de la gráfica (Opcinal, está por defecto en skyblue). (str)
  mu: Es la abreviación de GM (Opcional, trabajamos normalmente con mu=1). (float)
  """
  plt.plot(r, fa.ct2(r, C, True, mu=mu), color = col, linestyle= "--", linewidth=linew)
  plt.plot(r, fa.ct2(r, C, False, mu=mu), color = col, linestyle="-", linewidth=linew)
  return

def plotlightcone3(r, C, col="skyblue", mu = 1, zone2=True, linew=0.7):
  """
  Grafica las lineas de mundo de fotones entrantes y saleintes en coordenadas de Eddington-Filkelstein Retardadas.
  r: Distancias radiales que se desea graficar. (array)
  C: Constante arbitraría que desplaza hacia arriba o hacia abajo las curvas. (float)
  col: Color de la gráfica (Opcinal, está por defecto en skyblue). (str)
  mu: Es la abreviación de GM (Opcional, trabajamos normalmente con mu=1). (float)
  """
  plt.plot(r, fa.ct3(r, C, True, mu=mu), color = col, linestyle= "-", linewidth=linew)
  plt.plot(r, fa.ct3(r, C, False, mu=mu), color = col, linestyle="--", linewidth=linew)
  return


# Funciones de la sección 11.9 | Figura 11.6

ancho = 0.7

def hip_rconst(r, a, C, mu=1):
  col = "teal"
  if a >= 2*mu:
      plt.plot(fa.u_rphoton_mayor2mu_r(r, a, C, True), fa.v_rphoton_mayor2mu_r(r, a, C, True), color = col, linestyle="--", linewidth=ancho)
      plt.plot(-fa.u_rphoton_mayor2mu_r(r, a, C, True), fa.v_rphoton_mayor2mu_r(r, a, C, True), color = col, linestyle="--", linewidth=ancho)
  if a < 2*mu:
      plt.plot(fa.u_rphoton_menor2mu_r(r, a, C, True), fa.v_rphoton_menor2mu_r(r, a, C, True), color = col, linestyle="--", linewidth=ancho)
      plt.plot(fa.u_rphoton_menor2mu_r(r, a, C, True), -fa.v_rphoton_menor2mu_r(r, a, C, True), color = col, linestyle="--", linewidth=ancho)

      plt.plot(fa.u_rphoton_menor2mu_r(r, a, C, False), fa.v_rphoton_menor2mu_r(r, a, C, False), color = col, linestyle="--", linewidth=ancho)
      plt.plot(fa.u_rphoton_menor2mu_r(r, a, C, False), -fa.v_rphoton_menor2mu_r(r, a, C, False), color = col, linestyle="--", linewidth=ancho)
  return

def rec_tconst(r, t, C, mu = 1):
  col = "olive"
  plt.plot(fa.u_rphoton_menor2mu_t(r, t, C, True), fa.v_rphoton_menor2mu_t(r, t, C, True), color = col, linestyle="--", linewidth=ancho)
  plt.plot(fa.u_rphoton_menor2mu_t(r, t, C, True), -fa.v_rphoton_menor2mu_t(r, t, C, True), color = col, linestyle="--", linewidth=ancho)

  plt.plot(fa.u_rphoton_mayor2mu_t(r, t, C, True), fa.v_rphoton_mayor2mu_t(r, t, C, True), color = col, linestyle="--", linewidth=ancho)
  plt.plot(fa.u_rphoton_mayor2mu_t(r, t, C, True), -fa.v_rphoton_mayor2mu_t(r, t, C, True), color = col, linestyle="--", linewidth=ancho)

  plt.plot(-fa.u_rphoton_menor2mu_t(r, t, C, False), fa.v_rphoton_menor2mu_t(r, t, C, False), color = col, linestyle="--", linewidth=ancho)
  plt.plot(-fa.u_rphoton_menor2mu_t(r, t, C, False), -fa.v_rphoton_menor2mu_t(r, t, C, False), color = col, linestyle="--", linewidth=ancho)

  plt.plot(-fa.u_rphoton_mayor2mu_t(r, t, C, False), fa.v_rphoton_mayor2mu_t(r, t, C, False), color = col, linestyle="--", linewidth=ancho)
  plt.plot(-fa.u_rphoton_mayor2mu_t(r, t, C, False), -fa.v_rphoton_mayor2mu_t(r, t, C, False), color = col, linestyle="--", linewidth=ancho)
  return

def plot_masspart(r0, mu=1):
  r_1 = np.linspace(0,r0*mu,100)
  plt.plot(fa.u_rpart_menor2mu_r(r0, r_1, mu)[0], fa.v_rpart_menor2mu_r(r0, r_1, mu)[0], color = "salmon", linewidth=ancho)
  plt.plot(fa.u_rpart_mayor2mu_r(r0, r_1, mu)[0], fa.v_rpart_mayor2mu_r(r0, r_1, mu)[0], color = "salmon", linewidth=ancho)
  return

def plot_photon(r0, mu=1):
  r_1 = np.linspace(0,r0*mu,100)
  plt.plot(fa.u_rpart_menor2mu_ra(r0, r_1, mu)[0], fa.v_rpart_menor2mu_ra(r0, r_1, mu)[0], color = "goldenrod", linewidth=ancho)
  plt.plot(fa.u_rpart_mayor2mu_ra(r0, r_1, mu)[0], fa.v_rpart_mayor2mu_ra(r0, r_1, mu)[0], color = "goldenrod", linewidth=ancho)
  return

def v_coord(u): # Hiperbolas
  return np.sqrt(u**2 +1), -np.sqrt(u**2 +1)

def v_coord_1(u, C): # Rectas
  return u+C, -u+C