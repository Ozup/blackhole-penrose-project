import numpy as np
import matplotlib.pyplot as plt
from . import functions_blackholes as fbh


import importlib
importlib.reload(fbh)
c=1

# Funciones auxiliares que no se usan directamente en los notebooks, se usan en functions_blackholes.py

# Funciones de la sección 11.3 | Figura 11.1

def ct(r, C, signo, mu=1): # | Ecuación: NA | Capítulo:11 | Seccción:11.3 | Página 251 |
  """
  Ecuaciones de las líneas de mundo de fotones entrantes y salientes en coordenadas de Schwarzschild.
  Retorna un array de tiempo correspondiente a los r elegidos.
  r: Distancias radiales que se desea graficar. (array)
  C: Constante arbitraría que desplaza hacia arriba o hacia abajo las curvas. (float)
  signo: Signo para diferenciar si es saliente o entrante lo que se desea graficar. Saliente: True | Entrante: False. (1 or -1)
  mu: Es la abreviación de GM (Opcional, trabajamos normalmente con mu=1). (float)
  """
  if signo == True: n = 1 # Fotón Saliente (Outgoing)
  if signo == False: n = -1 # Fotón Entrante (Incoming)
  return n * (r + 2*mu*np.log(abs(r/(2*mu)-1))) + C

# Funciones de la sección 11.4 | Figura 11.2

def time_tao(r0, tao, mu):
  """
  Dado un tao de interés genera el r y t correspondientes.
  Retorna r, t, y mu.
  r0: Distancia radial donde t=0 y tao=0. (float)
  tao: Tiempo propio del que se desea hallar el r y t correspondiente. (float)
  mu: Es la abreviación de GM (Opcional, trabajamos normalmente con mu=1). (float)
  """
  c = 1
  r = (2*mu*c**2 * ( ((r0**3)/(2*mu*c**2))**0.5 - 3/2*tao )**2 )**(1/3)
  return r, fbh.time_11_4(r0, r, mu), mu

# Funciones de la sección 11.5 | Figura 11.3

def ct2(r, C, signo, mu=1):
  """
  Ecuaciones de las líneas de mundo de fotones entrantes y salientes en coordenadas de Eddington-Filkelstein Avanzadas.
  Retorna un array de tiempo correspondiente a los r elegidos.
  r: Distancias radiales que se desea graficar. (array)
  C: Constante arbitraría que desplaza hacia arriba o hacia abajo las curvas. (float)
  signo: Signo para diferenciar si es saliente o entrante lo que se desea graficar. Saliente: 1 | Entrante: -1. (1 or -1)
  mu: Es la abreviación de GM (Opcional, trabajamos normalmente con mu=1). (float)
  """
  if signo == True:
    ct =  -r + C
  if signo == False:
    ct = r + 4*mu*np.log(abs(r/(2*mu)-1)) + C
  return ct


def ct3(r, C, signo, mu=1):
  """
  Ecuaciones de las líneas de mundo de fotones entrantes y salientes en coordenadas de Eddington-Filkelstein Retardadas.
  Retorna un array de tiempo correspondiente a los r elegidos.
  r: Distancias radiales que se desea graficar. (array)
  C: Constante arbitraría que desplaza hacia arriba o hacia abajo las curvas. (float)
  signo: Signo para diferenciar si es saliente o entrante lo que se desea graficar. Saliente: 1 | Entrante: -1. (1 or -1)
  mu: Es la abreviación de GM (Opcional, trabajamos normalmente con mu=1). (float)
  """
  if signo == True:
    ct =  r + C
  if signo == False:
    ct = - r - 4*mu*np.log(abs(r/(2*mu)-1)) + C
  return ct

# Funciones de la sección 11.9 | Figura 11.6

###################################################
##       Funciones de u y v para r fijo          ##
###################################################

# Cuando r > 2mu
def v_rphoton_mayor2mu_r(r, a, C, signo, mu=1):
  """
  Retorna la coordenada v de un fotón, para un r = a.
  Válida para la zona: r>2mu.
  r: Distancias radiales que se desea graficar. (array)
  a: Valor de r constante que se desea graficar (float)
  C: Constante arbitraría que desplaza hacia arriba o hacia abajo las curvas. (float)
  signo: Signo para diferenciar si es saliente o entrante lo que se desea graficar. Saliente: 1 | Entrante: -1. (1 or -1)
  mu: Es la abreviación de GM (Opcional, trabajamos normalmente con mu=1). (float)
  (Esta función grafica hiperbolas)
  """
  r = np.linspace(2*mu, r[-1],10000)
  t = ct(r, C, signo, mu)
  r = a
  return (r/(2*mu)-1)**(1/2) * np.exp(r/(4*mu)) * np.sinh(c*t/(4*mu))

def u_rphoton_mayor2mu_r(r, a, C, signo, mu=1):
  """
  Retorna la coordenada u de un fotón, para un r = a.
  Válida para la zona: r>2mu.
  r: Distancias radiales que se desea graficar. (array)
  a: Valor de r constante que se desea graficar (float)
  C: Constante arbitraría que desplaza hacia arriba o hacia abajo las curvas. (float)
  signo: Signo para diferenciar si es saliente o entrante lo que se desea graficar. Saliente: 1 | Entrante: -1. (1 or -1)
  mu: Es la abreviación de GM (Opcional, trabajamos normalmente con mu=1). (float)
  (Esta función grafica hiperbolas)
  """
  r = np.linspace(2*mu, r[-1],10000)
  t = ct(r, C, signo, mu)
  r = a
  return (r/(2*mu)-1)**(1/2) * np.exp(r/(4*mu)) * np.cosh(c*t/(4*mu))

# Cuando r < 2mu
def v_rphoton_menor2mu_r(r, a, C, signo, mu=1):
  """
  Retorna la coordenada v de un fotón, para un r = a.
  Válida para la zona: r<2mu.
  r: Distancias radiales que se desea graficar. (array)
  a: Valor de r constante que se desea graficar (float)
  C: Constante arbitraría que desplaza hacia arriba o hacia abajo las curvas. (float)
  signo: Signo para diferenciar si es saliente o entrante lo que se desea graficar. Saliente: 1 | Entrante: -1. (1 or -1)
  mu: Es la abreviación de GM (Opcional, trabajamos normalmente con mu=1). (float)
  (Esta función grafica hiperbolas)
  """
  r = np.linspace(r[0], 2*mu,10000)
  t = ct(r, C, signo, mu)
  r = a
  return (1-r/(2*mu))**(1/2) * np.exp(r/(4*mu)) * np.cosh(c*t/(4*mu))

def u_rphoton_menor2mu_r(r, a, C, signo, mu=1):
  """
  Retorna la coordenada u de un fotón, para un r = a.
  Válida para la zona: r<2mu.
  r: Distancias radiales que se desea graficar. (array)
  a: Valor de r constante que se desea graficar (float)
  C: Constante arbitraría que desplaza hacia arriba o hacia abajo las curvas. (float)
  signo: Signo para diferenciar si es saliente o entrante lo que se desea graficar. Saliente: 1 | Entrante: -1. (1 or -1)
  mu: Es la abreviación de GM (Opcional, trabajamos normalmente con mu=1). (float)
  (Esta función grafica hiperbolas)
  """
  r = np.linspace(r[0], 2*mu,10000)
  t = ct(r, C, signo, mu)
  r = a

  return (1-r/(2*mu))**(1/2) * np.exp(r/(4*mu)) * np.sinh(c*t/(4*mu))

###################################################
##       Funciones de u y v para t fijo          ##
###################################################

# Cuando r > 2mu
def v_rphoton_mayor2mu_t(r, t, C, signo, mu=1):
  """
  Retorna la coordenada v de un fotón, para un t = t.
  Válida para la zona: r>2mu.
  r: Distancias radiales que los t que se desea graficar. (array)
  t: Valor de t constante que se desea graficar (float)
  C: Constante arbitraría que desplaza hacia arriba o hacia abajo las curvas. (float)
  signo: Signo para diferenciar si es saliente o entrante lo que se desea graficar. Saliente: 1 | Entrante: -1. (1 or -1)
  mu: Es la abreviación de GM (Opcional, trabajamos normalmente con mu=1). (float)
  (Esta función grafica rectas)
  """
  r = np.linspace(2*mu, r[-1],10000)
  return (r/(2*mu)-1)**(1/2) * np.exp(r/(4*mu)) * np.sinh(c*t/(4*mu))

def u_rphoton_mayor2mu_t(r, t, C, signo, mu=1):
  """
  Retorna la coordenada u de un fotón, para un t = t.
  Válida para la zona: r>2mu.
  r: Distancias radiales que los t que se desea graficar. (array)
  t: Valor de t constante que se desea graficar (float)
  C: Constante arbitraría que desplaza hacia arriba o hacia abajo las curvas. (float)
  signo: Signo para diferenciar si es saliente o entrante lo que se desea graficar. Saliente: 1 | Entrante: -1. (1 or -1)
  mu: Es la abreviación de GM (Opcional, trabajamos normalmente con mu=1). (float)
  (Esta función grafica rectas)
  """
  r = np.linspace(2*mu, r[-1],10000)
  return (r/(2*mu)-1)**(1/2) * np.exp(r/(4*mu)) * np.cosh(c*t/(4*mu))

# Cuando r < 2mu

def v_rphoton_menor2mu_t(r, t, C, signo, mu=1):
  """
  Retorna la coordenada v de un fotón, para un t = t.
  Válida para la zona: r<2mu.
  r: Distancias radiales que los t que se desea graficar. (array)
  t: Valor de t constante que se desea graficar (float)
  C: Constante arbitraría que desplaza hacia arriba o hacia abajo las curvas. (float)
  signo: Signo para diferenciar si es saliente o entrante lo que se desea graficar. Saliente: 1 | Entrante: -1. (1 or -1)
  mu: Es la abreviación de GM (Opcional, trabajamos normalmente con mu=1). (float)
  (Esta función grafica rectas)
  """
  r = np.linspace(r[0], 2*mu,10000)
  return (1-r/(2*mu))**(1/2) * np.exp(r/(4*mu)) * np.cosh(c*t/(4*mu))

def u_rphoton_menor2mu_t(r, t, C, signo, mu=1):
  """
  Retorna la coordenada u de un fotón, para un t = t.
  Válida para la zona: r<2mu.
  r: Distancias radiales que los t que se desea graficar. (array)
  t: Valor de t constante que se desea graficar (float)
  C: Constante arbitraría que desplaza hacia arriba o hacia abajo las curvas. (float)
  signo: Signo para diferenciar si es saliente o entrante lo que se desea graficar. Saliente: 1 | Entrante: -1. (1 or -1)
  mu: Es la abreviación de GM (Opcional, trabajamos normalmente con mu=1). (float)
  (Esta función grafica rectas)
  """
  r = np.linspace(r[0], 2*mu,10000)
  return (1-r/(2*mu))**(1/2) * np.exp(r/(4*mu)) * np.sinh(c*t/(4*mu))

###################################################
###################################################
##         -= PARTÍCULAS CON MASA =-             ##
###################################################
###################################################

# Cuando r > 2mu
def v_rpart_mayor2mu_r(r0, r, mu):
  r = np.linspace(2*mu, r[-1],1000)[::-1]
  t = fbh.time_11_4(r0, r, mu=1)
  v = (r/(2*mu)-1)**(1/2) * np.exp(r/(4*mu)) * np.sinh(c*t/(4*mu))
  return v, r, t

def u_rpart_mayor2mu_r(r0, r, mu):
  r = np.linspace(2*mu, r[-1],1000)[::-1]
  t = fbh.time_11_4(r0, r, mu=1)
  u = (r/(2*mu)-1)**(1/2) * np.exp(r/(4*mu)) * np.cosh(c*t/(4*mu))
  return u, r, t

# Cuando r < 2mu
def v_rpart_menor2mu_r(r0, r, mu):
  r = np.linspace(r[0], 2*mu,1000)[::-1]
  t = fbh.time_11_4(r0, r, mu=1)
  v = (1-r/(2*mu))**(1/2) * np.exp(r/(4*mu)) * np.cosh(c*t/(4*mu))
  return v, r, t

def u_rpart_menor2mu_r(r0, r, mu):
  r = np.linspace(r[0], 2*mu,1000)[::-1]
  t = fbh.time_11_4(r0, r, mu=1)
  u = (1-r/(2*mu))**(1/2) * np.exp(r/(4*mu)) * np.sinh(c*t/(4*mu))
  return u, r, t

###################################################
###################################################
##                -= FOTONES =-                 ##
###################################################
###################################################

# Cuando r > 2mu
def v_rpart_mayor2mu_ra(r0, r, a, mu=1):
  r = np.linspace(2*mu, r[-1],1000)[::-1]
  signo = False
  t = -r-2*np.log(abs(r/2-1))
  v = (r/2-1)**(1/2) * np.exp(r/4) * np.sinh(t/4)
  return v, r, t

def u_rpart_mayor2mu_ra(r0, r, a, mu=1):
  r = np.linspace(2*mu, r[-1],1000)[::-1]
  signo = False
  t = -r-2*np.log(abs(r/2-1))
  u = (r/2-1)**(1/2) * np.exp(r/4) * np.cosh(t/4)
  return u, r, t

# Cuando r < 2mu
def v_rpart_menor2mu_ra(r0, r, a, mu=1):
  r = np.linspace(r[0], 2*mu,1000)[::-1]
  signo = False
  t = -r-2*np.log(abs(r/2-1))
  v = (1-r/2)**(1/2) * np.exp(r/4) * np.cosh(t/4)
  return v, r, t

def u_rpart_menor2mu_ra(r0, r, a, mu=1):
  r = np.linspace(r[0], 2*mu,1000)[::-1]
  signo = False
  t = -r-2*np.log(abs(r/2-1))
  u = (1-r/2)**(1/2) * np.exp(r/4) * np.sinh(t/4)
  return u, r, t

# ----------------------------------------------#