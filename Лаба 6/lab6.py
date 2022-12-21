import numpy as np 
import matplotlib.pyplot as plt
from slae import *

def rabbit(a, x):
  return a*x*(1 - x)

def diff_r(a, x):
  return a - 2*a*x

def for_N(a, x, tau, xk):
  return x - a*tau*x*(1-x) - xk

def diff_N(a, x, tau):
  return 1 - tau*a*(1-x) + tau*a*x

def NewtonEq(n, x, a, tau, xk):
  for i in range(n):
    x = x - for_N(a, x, tau, xk)/diff_N(a, x, tau)
  return x

def for_Nt(a, x, tau, xk):
  return x - tau*(rabbit(a, x) + rabbit(a, xk))/2 - xk

def diff_Nt(a, x, tau, xk):
  return 1 - tau*(diff_r(a, x) + diff_r(a, xk))

def NewtonEqt(n, x, a, tau, xk):
  for i in range(n):
    x = x - for_Nt(a, x, tau, xk)/diff_Nt(a, x, tau, xk)
  return x

def ans_r(a, t):
  return np.exp(a*t)/(np.exp(a*t) + 10)

def euc(x, y):
  return np.sqrt(pow(x[0] - y[0], 2) + pow(x[1] - y[1], 2))

#========= Явная схема Эйлера уравнение ===============

Err = list()
grid = list()
a = 5
cond = ans_r(a, 0)
for i in range(1, 15):
  Tsls = np.linspace(0, 2, pow(2, i))
  xj = np.array(cond)
  tau = Tsls[1]
  eps = 0
  for j in range(1, len(Tsls)):
    xj = xj + rabbit(a, xj)*tau
    e = float(np.abs(ans_r(a, Tsls[j]) - xj))
    if e > eps:
      eps = e
  Err.append(eps)
  grid.append(tau)

print("Явная схема Эйлера для уравнения")

plt.figure(1, figsize=(10, 5))
plt.subplot(121)
plt.autoscale(tight=True)
plt.ylabel(r'$y$')
plt.xlabel(r'$\tau$')
plt.title("Явная схема Эйлера для уравнения")
Tsls = np.linspace(0, 2, pow(2, 4))
xj = np.array(cond)
tau = Tsls[1]
y = [xj]
for j in range(1, len(Tsls)):
  xj = xj + rabbit(a, xj)*tau
  y.append(xj)
t = np.linspace(0, 2, 1000)
plt.plot(t, ans_r(a, t), linewidth=1, color='blue')
plt.plot(Tsls, y, linewidth=1, color='orange')
plt.grid()
plt.subplot(122)
plt.autoscale(tight=True)
plt.ylabel(r'$\delta$')
plt.xlabel(r'$\tau$')
plt.title("Явная схема Эйлера для уравнения")
E_x = np.log(grid)
E_y = np.log(Err)
plt.plot(E_x, E_y, linewidth=1, color='blue')
plt.errorbar(E_x, E_y, yerr=0.00, xerr=0.00, fmt='.', ecolor='red', color='red')
plt.grid()
print("Порядок сходимости = ", (E_y[2] - E_y[-1])/(E_x[2] - E_x[-1]))
print()

#========= Неявная схема Эйлера ===============

Err = list()
grid = list()
a = 5
cond = ans_r(a, 0)
for i in range(1, 10):
  Tsls = np.linspace(0, 2, pow(2, i))
  xj = np.array(cond)
  tau = Tsls[1]
  eps = 0
  for j in range(1, len(Tsls)):
    xj = NewtonEq(50, xj, a, tau, xj)
    e = float(np.abs(ans_r(a, Tsls[j]) - xj))
    if e > eps:
      eps = e
  Err.append(eps)
  grid.append(tau)

print("Неявная схема Эйлера")

plt.figure(2, figsize=(10, 5))
plt.subplot(121)
plt.autoscale(tight=True)
plt.ylabel(r'$y$')
plt.xlabel(r'$\tau$')
plt.title("Неявная схема Эйлера для уравнения")
Tsls = np.linspace(0, 2, pow(2, 4))
xj = np.array(cond)
tau = Tsls[1]
y = [xj]
for j in range(1, len(Tsls)):
  xj = NewtonEq(50, xj, a, tau, xj)
  y.append(xj)
t = np.linspace(0, 2, 1000)
plt.plot(t, ans_r(a, t), linewidth=1, color='blue')
plt.plot(Tsls, y, linewidth=1, color='orange')
plt.grid()
plt.subplot(122)
plt.autoscale(tight=True)
plt.ylabel(r'$\delta$')
plt.xlabel(r'$\tau$')
plt.title("Неявная схема Эйлера для уравнения")
E_x = np.log(grid)
E_y = np.log(Err)
plt.plot(E_x, E_y, linewidth=1, color='blue')
plt.errorbar(E_x, E_y, yerr=0.00, xerr=0.00, fmt='.', ecolor='red', color='red')
plt.grid()
print("Порядок сходимости = ", (E_y[4] - E_y[-1])/(E_x[4] - E_x[-1]))
print()

#========= Метод трапеций ===============

Err = list()
grid = list()
a = 5
cond = ans_r(a, 0)
for i in range(1, 15):
  Tsls = np.linspace(0, 2, pow(2, i))
  xj = np.array(cond)
  tau = Tsls[1]
  eps = 0
  for j in range(1, len(Tsls)):
    xj = NewtonEqt(50, xj, a, tau, xj)
    e = float(np.abs(ans_r(a, Tsls[j]) - xj))
    if e > eps:
      eps = e
  Err.append(eps)
  grid.append(tau)

print("Метод трапеций")

plt.figure(3, figsize=(10, 5))
plt.subplot(121)
plt.autoscale(tight=True)
plt.ylabel(r'$y$')
plt.xlabel(r'$\tau$')
plt.title("Метод трапеций")
Tsls = np.linspace(0, 2, pow(2, 5))
xj = np.array(cond)
tau = Tsls[1]
y = [xj]
for j in range(1, len(Tsls)):
  xj = NewtonEqt(50, xj, a, tau, xj)
  y.append(xj)
t = np.linspace(0, 2, 1000)
plt.plot(t, ans_r(a, t), linewidth=1, color='blue')
plt.plot(Tsls, y, linewidth=1, color='orange')
plt.grid()
plt.subplot(122)
plt.autoscale(tight=True)
plt.ylabel(r'$\delta$')
plt.xlabel(r'$\tau$')
plt.title("Метод трапеций")
E_x = np.log(grid)
E_y = np.log(Err)
plt.plot(E_x, E_y, linewidth=1, color='blue')
plt.errorbar(E_x, E_y, yerr=0.00, xerr=0.00, fmt='.', ecolor='red', color='red')
plt.grid()
print("Порядок сходимости = ", (E_y[5] - E_y[-1])/(E_x[5] - E_x[-1]))
print()

#========= Метод РК2 ===============

Err = list()
grid = list()
a = 5
cond = ans_r(a, 0)
for i in range(1, 15):
  Tsls = np.linspace(0, 2, pow(2, i))
  xj = np.array(cond)
  tau = Tsls[1]
  eps = 0
  for j in range(1, len(Tsls)):
    xj = xj + tau*(rabbit(a, xj) + rabbit(a, xj + tau*rabbit(a, xj)))/2
    e = float(np.abs(ans_r(a, Tsls[j]) - xj))
    if e > eps:
      eps = e
  Err.append(eps)
  grid.append(tau)

print("Метод РК2")

plt.figure(4, figsize=(10, 5))
plt.subplot(121)
plt.autoscale(tight=True)
plt.ylabel(r'$y$')
plt.xlabel(r'$\tau$')
plt.title("Метод РК2")
Tsls = np.linspace(0, 2, pow(2, 5))
xj = np.array(cond)
tau = Tsls[1]
y = [xj]
for j in range(1, len(Tsls)):
  xj = xj + tau*(rabbit(a, xj) + rabbit(a, xj + tau*rabbit(a, xj)))/2
  y.append(xj)
t = np.linspace(0, 2, 1000)
plt.plot(t, ans_r(a, t), linewidth=1, color='blue')
plt.plot(Tsls, y, linewidth=1, color='orange')
plt.grid()
plt.subplot(122)
plt.autoscale(tight=True)
plt.ylabel(r'$\delta$')
plt.xlabel(r'$\tau$')
plt.title("Метод РК2")
E_x = np.log(grid)
E_y = np.log(Err)
plt.plot(E_x, E_y, linewidth=1, color='blue')
plt.errorbar(E_x, E_y, yerr=0.00, xerr=0.00, fmt='.', ecolor='red', color='red')
plt.grid()
print("Порядок сходимости = ", (E_y[5] - E_y[-1])/(E_x[5] - E_x[-1]))
print()

#========= Метод РК4 ===============

Err = list()
grid = list()
a = 5
cond = ans_r(a, 0)
for i in range(1, 15):
  Tsls = np.linspace(0, 2, pow(2, i))
  xj = np.array(cond)
  tau = Tsls[1]
  eps = 0
  for j in range(1, len(Tsls)):
    k1 = rabbit(a, xj)
    k2 = rabbit(a, xj + tau/2*k1)
    k3 = rabbit(a, xj + tau/2*k2)
    k4 = rabbit(a, xj + tau*k3)
    xj = xj + tau*(k1+2*k2+2*k3+k4)/6
    e = float(np.abs(ans_r(a, Tsls[j]) - xj))
    if e > eps:
      eps = e
  Err.append(eps)
  grid.append(tau)

print("Метод РК4")

plt.figure(5, figsize=(10, 5))
plt.subplot(121)
plt.autoscale(tight=True)
plt.ylabel(r'$y$')
plt.xlabel(r'$\tau$')
plt.title("Метод РК4")
Tsls = np.linspace(0, 2, pow(2, 5))
xj = np.array(cond)
tau = Tsls[1]
y = [xj]
for j in range(1, len(Tsls)):
  k1 = rabbit(a, xj)
  k2 = rabbit(a, xj + tau/2*k1)
  k3 = rabbit(a, xj + tau/2*k2)
  k4 = rabbit(a, xj + tau*k3)
  xj = xj + tau*(k1+2*k2+2*k3+k4)/6
  y.append(xj)
t = np.linspace(0, 2, 1000)
plt.plot(t, ans_r(a, t), linewidth=1, color='blue')
plt.plot(Tsls, y, linewidth=1, color='orange')
plt.grid()
plt.subplot(122)
plt.autoscale(tight=True)
plt.ylabel(r'$\delta$')
plt.xlabel(r'$\tau$')
plt.title("Метод РК4")
E_x = np.log(grid)
E_y = np.log(Err)
plt.plot(E_x, E_y, linewidth=1, color='blue')
plt.errorbar(E_x, E_y, yerr=0.00, xerr=0.00, fmt='.', ecolor='red', color='red')
plt.grid()
print("Порядок сходимости = ", (E_y[3] - E_y[-3])/(E_x[3] - E_x[-3]))
print()

#============== Модель Лотки-Вольтерры =================

def Lotki(a, b, c, d, x):
  res = np.zeros((2,1))
  res[0] = (a - b*x[1])*x[0]
  res[1] = (-c + d*x[0])*x[1]
  return res

a = 0.5
b = 0.2
c = 0.1
d = 0.2

plt.figure(6, figsize=(10, 5))
plt.subplot(121)
plt.autoscale(tight=True)
plt.ylabel(r'$N_x$')
plt.xlabel(r'$N_y$')
plt.title("Метод РК4")
Tsls = np.linspace(0, 100, 10*pow(2, 5))
tau = Tsls[1]
xj = np.ones((2,1))
xj[0] = 0.95
xj[1] = 0.8
x=[float(xj[0])]
y=[float(xj[1])]
for j in range(1, len(Tsls)):
  k1 = Lotki(a, b, c, d, xj)
  k2 = Lotki(a, b, c, d, xj + tau/2*k1)
  k3 = Lotki(a, b, c, d, xj + tau/2*k2)
  k4 = Lotki(a, b, c, d, xj + tau*k3)
  xj = xj + tau*(k1+2*k2+2*k3+k4)/6
  x.append(float(xj[0]))
  y.append(float(xj[1]))
plt.plot(x, y, linewidth=1, color='blue')
plt.grid()

plt.subplot(122)
plt.autoscale(tight=True)
plt.ylabel(r'$N$')
plt.xlabel(r'$\tau$')
plt.title("Модель Лотки-Вольтеры")
plt.plot(Tsls, x, linewidth=1, color='blue')
plt.plot(Tsls, y, linewidth=1, color='red')
plt.grid()

#================ Система Лоренца =================

def Lorenz(s, r, b, x):
  res = np.zeros((3,1))
  res[0] = s*(x[1] - x[0])
  res[1] = x[0]*(r - x[2]) - x[1]
  res[2] = x[0]*x[1] - b*x[2]
  return res

s = 10
r = 28
b = 8.0/3

Tsls = np.linspace(0, 50, 100*pow(2, 5))
tau = Tsls[1]
xj = np.ones((3,1))
v = np.zeros((3,1))
v[0] = 1
eps = []
xj[0] = 3.051522
xj[1] = 1.582542
xj[2] = 15.62388
x=[float(xj[0])]
y=[float(xj[1])]
z=[float(xj[2])]
for j in range(1, len(Tsls)):
  k1 = Lorenz(s, r, b, xj)
  k2 = Lorenz(s, r, b, xj + tau/2*k1)
  k3 = Lorenz(s, r, b, xj + tau/2*k2)
  k4 = Lorenz(s, r, b, xj + tau*k3)
  xj = xj + tau*(k1+2*k2+2*k3+k4)/6
  k1 = Lorenz(s, r, b, v)
  k2 = Lorenz(s, r, b, v + tau/2*k1)
  k3 = Lorenz(s, r, b, v + tau/2*k2)
  k4 = Lorenz(s, r, b, v + tau*k3)
  v = v + tau*(k1+2*k2+2*k3+k4)/6
  eps.append(float(np.sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])))
  v = v/eps[-1]
  x.append(float(xj[0]))
  y.append(float(xj[1]))
  z.append(float(xj[2]))
eps = np.array(eps)
print("Показатель Ляпунова = ", sum(np.log(eps))/50)

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
ax.plot(x, y, z, lw=0.5)
ax.set_xlabel("X Axis")
ax.set_ylabel("Y Axis")
ax.set_zlabel("Z Axis")
ax.set_title("Система Лоренца, r=28")

#============= Бифуркации =====================

s = 10
r = 1
b = 8.0/3

Tsls = np.linspace(0, 50, 100*pow(2, 5))
tau = Tsls[1]
xj = np.ones((3,1))


plt.figure(figsize=(5, 5))
plt.autoscale(tight=True)
plt.ylabel(r'$intersections$')
plt.xlabel(r'$r$')
plt.title("Бифуркационная диаграмма")
plt.grid()
rlist = np.linspace(1, 30, 90)
z0 = 12
for r in rlist:
  bif = []
  xj[0] = 3.051522
  xj[1] = 1.582542
  xj[2] = 15.62388
  x = float(xj[0])
  y = float(xj[1])
  z = float(xj[2])
  for j in range(1, len(Tsls)):
    k1 = Lorenz(s, r, b, xj)
    k2 = Lorenz(s, r, b, xj + tau/2*k1)
    k3 = Lorenz(s, r, b, xj + tau/2*k2)
    k4 = Lorenz(s, r, b, xj + tau*k3)
    xj = xj + tau*(k1+2*k2+2*k3+k4)/6
    if (z < z0 and xj[2] > z0) or (z > z0 and xj[2] < z0):
      bif.append(float(x + np.abs((z-z0)/(z-xj[2]))*(xj[0] - x)))
    x = float(xj[0])
    y = float(xj[1])
    z = float(xj[2])
  plt.errorbar([r for i in bif], bif, yerr=0.00, xerr=0.00, fmt='.', ecolor='red', color='red')


plt.show()