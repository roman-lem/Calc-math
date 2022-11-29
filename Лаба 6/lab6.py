import numpy as np 
import matplotlib.pyplot as plt
from slae import *

def func(x):
  res = np.zeros((2, 1))
  res[0, 0] = x[1] - 2*x[1]
  res[1, 0] = 5*x[0] - x[1]
  return res

def ans(t):
  res = np.zeros((2, 1))
  res[0, 0] = np.cos(3*t)
  res[1, 0] = (np.cos(3*t) + 3*np.sin(3*t))/2
  return res

def Atau(tau):
  ret = np.zeros((2, 2))
  ret[0, 0] = 1 - tau
  ret[0, 1] = 2 * tau
  ret[1, 0] = -5 * tau
  ret[1, 1] = 1 + tau
  return ret

def euc(x, y):
  return np.sqrt(pow(x[0] - y[0], 2) + pow(x[1] - y[1], 2))

Cond = ans(0)

#========= Явная схема Эйлера ===============

Err = list()
grid = list()
for i in range(1, 15):
  Tsls = np.linspace(0, 1, pow(2, i))
  xj = np.array(Cond)
  tau = Tsls[1]
  eps = 0
  for j in range(1, len(Tsls)):
    xj = xj + func(xj)*tau
    e = int(euc(ans(Tsls[j]), xj))
    if e > eps:
      eps = e
  Err.append(eps)
  grid.append(tau)

print("Явная схема Эйлера")
print(ans(1))
print("delta = ", Err[-1])
print()

plt.figure(1, figsize=(10, 5))
plt.subplot(121)
plt.autoscale(tight=True)
plt.ylabel(r'$\delta$')
plt.xlabel(r'$\tau$')
plt.title("Явная схема Эйлера")
E_x = (grid)
E_y = (Err)
plt.plot(E_x, E_y, linewidth=1, color='blue')
plt.errorbar(E_x, E_y, yerr=0.00, xerr=0.00, fmt='.', ecolor='red', color='red')
plt.grid()

#========= Неявная схема Эйлера ===============

Err = list()
grid = list()
for i in range(1, 15):
  Tsls = np.linspace(0, 1, pow(2, i))
  xj = np.array(Cond)
  tau = Tsls[1]
  eps = 0
  for j in range(1, len(Tsls)):
    xj = Hauss(Atau(tau), xj, 2)
    e = int(euc(ans(Tsls[j]), xj))
    if e > eps:
      eps = e
  Err.append(eps)
  grid.append(tau)

print("Неявная схема Эйлера")
print(ans(1))
print("delta = ", Err[-1])
print()

plt.subplot(122)
plt.autoscale(tight=True)
plt.ylabel(r'$\delta$')
plt.xlabel(r'$\tau$')
plt.title("Неявная схема Эйлера")
E_x = (grid)
E_y = (Err)
plt.plot(E_x, E_y, linewidth=1, color='blue')
plt.errorbar(E_x, E_y, yerr=0.00, xerr=0.00, fmt='.', ecolor='red', color='red')
plt.grid()


plt.show()