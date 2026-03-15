import numpy as np
import matplotlib.pyplot as plt
from . import functions_blackholes as fbh
from . import functions_aux as fa


import importlib
importlib.reload(fbh)
c=1

## Penrose TF ##

# Rotation
def v(r, t):
  return t + r
def w(r, t):
  return t - r
  
# Compactification 
def p(r, t):
  return np.arctan(v(r,t))
def q(r, t):
  return np.arctan(w(r,t))

# Rotation
def tprime(r, t):
  return p(r, t) + q(r, t)
def rprime(r, t):
  return p(r, t) - q(r, t)

# Funciones para graficar Penrose en Kruskal

C = 0
def curv_rcons_menor(r, a, color = "black", wi = 0.5, line = "--"):
  A1 = fa.u_rphoton_menor2mu_r(r, a, C, True)
  B1 = fa.v_rphoton_menor2mu_r(r, a, C, True)

  plt.plot(rprime(A1,B1), tprime(A1,B1), color = color, linewidth=wi, linestyle=line)
  plt.plot(rprime(A1,-B1), tprime(A1,-B1), color = color, linewidth=wi, linestyle=line)

  plt.plot(rprime(-A1,B1), tprime(-A1,B1), color = color, linewidth=wi, linestyle=line)
  plt.plot(rprime(-A1,-B1), tprime(-A1,-B1), color = color, linewidth=wi, linestyle=line)
  return

def curv_rcons_mayor(r, a, color = "black", wi = 0.5, line = "--"):
  A2 = fa.u_rphoton_mayor2mu_r(r, a, C, True)
  B2 = fa.v_rphoton_mayor2mu_r(r, a, C, True)


  plt.plot(rprime(A2,B2), tprime(A2,B2), color = color, linewidth=wi, linestyle=line)
  plt.plot(rprime(A2,-B2), tprime(A2,-B2), color = color, linewidth=wi, linestyle=line)
  plt.plot(rprime(-A2,B2), tprime(-A2,B2), color = color, linewidth=wi, linestyle=line)
  plt.plot(rprime(-A2,-B2), tprime(-A2,-B2), color = color, linewidth=wi, linestyle=line)
  return

def curv_tcons(r, t, color = "black", wi = 0.5, line = "--"):
  A = fa.u_rphoton_menor2mu_t(r, t, C, True)
  B = fa.v_rphoton_menor2mu_t(r, t, C, True)

  plt.plot(rprime(A, B), tprime(A, B), color = color, linestyle=line, linewidth=wi)
  plt.plot(rprime(A, -B), tprime(A, -B), color = color, linestyle=line, linewidth=wi)

  A = fa.u_rphoton_mayor2mu_t(r, t, C, True)
  B = fa.v_rphoton_mayor2mu_t(r, t, C, True)
  plt.plot(rprime(A, B), tprime(A, B), color = color, linestyle=line, linewidth=wi)
  plt.plot(rprime(A, -B), tprime(A, -B), color = color, linestyle=line, linewidth=wi)

  A = fa.u_rphoton_menor2mu_t(r, t, C, False)
  B = fa.v_rphoton_menor2mu_t(r, t, C, False)
  plt.plot(rprime(A, B), tprime(A, B), color = color, linestyle=line, linewidth=wi) # Sobra aquí.
  plt.plot(rprime(-A, B), tprime(-A, B), color = color, linestyle=line, linewidth=wi)
  plt.plot(rprime(-A, -B), tprime(-A, -B), color = color, linestyle=line, linewidth=wi)


  A = fa.u_rphoton_mayor2mu_t(r, t, C, False)
  B = fa.v_rphoton_mayor2mu_t(r, t, C, False)
  plt.plot(rprime(-A, B), tprime(-A, B), color = color, linestyle=line, linewidth=wi)
  plt.plot(rprime(-A, -B), tprime(-A, -B), color = color, linestyle=line, linewidth=wi)
  return


def plot_masspart_penrose(r0, mu=1):
  r_1 = np.linspace(0,r0*mu,100)
  A = fa.u_rpart_menor2mu_r(r0, r_1, mu)[0]
  B = fa.v_rpart_menor2mu_r(r0, r_1, mu)[0]
  plt.plot(rprime(A, B), tprime(A, B), color = "black", linewidth=0.7)
  A = fa.u_rpart_mayor2mu_r(r0, r_1, mu)[0]
  B = fa.v_rpart_mayor2mu_r(r0, r_1, mu)[0]
  plt.plot(rprime(A, B), tprime(A, B), color = "black", linewidth=0.7)
  return

def plot_photon_penrose(r0, mu=1):


  r_1 = np.linspace(0,r0*mu,100)
  A = fa.u_rpart_menor2mu_ra(r0, r_1, mu)[0]
  B = fa.v_rpart_menor2mu_ra(r0, r_1, mu)[0]
  plt.plot(rprime(A, B), tprime(A, B), color = "black", linewidth=0.7)
  A = fa.u_rpart_mayor2mu_ra(r0, r_1, mu)[0]
  B = fa.v_rpart_mayor2mu_ra(r0, r_1, mu)[0]
  plt.plot(rprime(A, B), tprime(A, B), color = "black", linewidth=0.7)
  return



#Modificaciones para translación

def curv_rcons_menor_T(r, a, T, color = "black", wi = 0.5, line = "--"):
  A1 = fa.u_rphoton_menor2mu_r(r, a, C, True)
  B1 = fa.v_rphoton_menor2mu_r(r, a, C, True)

  plt.plot(rprime(A1,B1), tprime(A1,B1) + T, color = color, linewidth=wi, linestyle=line)
  plt.plot(rprime(A1,-B1), tprime(A1,-B1) + T, color = color, linewidth=wi, linestyle=line)

  plt.plot(rprime(-A1,B1), tprime(-A1,B1) + T, color = color, linewidth=wi, linestyle=line)
  plt.plot(rprime(-A1,-B1), tprime(-A1,-B1) + T, color = color, linewidth=wi, linestyle=line)
  return

def curv_rcons_mayor_T(r, a, T, color = "black", wi = 0.5, line = "--"):
  A2 = fa.u_rphoton_mayor2mu_r(r, a, C, True)
  B2 = fa.v_rphoton_mayor2mu_r(r, a, C, True)


  plt.plot(rprime(A2,B2), tprime(A2,B2) + T, color = color, linewidth=wi, linestyle=line)
  plt.plot(rprime(A2,-B2), tprime(A2,-B2) + T, color = color, linewidth=wi, linestyle=line)
  plt.plot(rprime(-A2,B2), tprime(-A2,B2) + T, color = color, linewidth=wi, linestyle=line)
  plt.plot(rprime(-A2,-B2), tprime(-A2,-B2) + T, color = color, linewidth=wi, linestyle=line)
  return

def curv_tcons_T(r, t, T, color = "black", wi = 0.5, line = "--"):
  A = fa.u_rphoton_menor2mu_t(r, t, C, True)
  B = fa.v_rphoton_menor2mu_t(r, t, C, True)

  plt.plot(rprime(A, B), tprime(A, B) + T, color = color, linestyle=line, linewidth=wi)
  plt.plot(rprime(A, -B), tprime(A, -B) + T, color = color, linestyle=line, linewidth=wi)

  A = fa.u_rphoton_mayor2mu_t(r, t, C, True)
  B = fa.v_rphoton_mayor2mu_t(r, t, C, True)
  plt.plot(rprime(A, B), tprime(A, B) + T, color = color, linestyle=line, linewidth=wi)
  plt.plot(rprime(A, -B), tprime(A, -B) + T, color = color, linestyle=line, linewidth=wi)

  A = fa.u_rphoton_menor2mu_t(r, t, C, False)
  B = fa.v_rphoton_menor2mu_t(r, t, C, False)
  plt.plot(rprime(A, B), tprime(A, B) + T, color = color, linestyle=line, linewidth=wi) # Sobra aquí.
  plt.plot(rprime(-A, B), tprime(-A, B) + T, color = color, linestyle=line, linewidth=wi)
  plt.plot(rprime(-A, -B), tprime(-A, -B) + T, color = color, linestyle=line, linewidth=wi)


  A = fa.u_rphoton_mayor2mu_t(r, t, C, False)
  B = fa.v_rphoton_mayor2mu_t(r, t, C, False)
  plt.plot(rprime(-A, B), tprime(-A, B) + T, color = color, linestyle=line, linewidth=wi)
  plt.plot(rprime(-A, -B), tprime(-A, -B) + T, color = color, linestyle=line, linewidth=wi)
  return

