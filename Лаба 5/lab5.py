import numpy as np 
import matplotlib.pyplot as plt
import random as rd

def func(x):
  return 5*pow(x, 4)

ans = 5**5

def rect(xk, fk):
  res = 0
  for i in range(len(xk)-1):
    res += fk[i]*(xk[i+1] - xk[i])
  return res

def trape(xk, fk):
  res = 0
  for i in range(len(xk)-1):
    res += (fk[i+1] + fk[i])/2*(xk[i+1] - xk[i])
  return res

def simp(xk, fk):
  res = 0
  for i in range(1, len(xk)-1, 2):
    res += (xk[i+1] - xk[i-1])/6*(fk[i-1] + 4*fk[i] + fk[i+1])
  return res

def rect_R(xk, fk):
  res = 0
  for i in range(len(xk)-1):
    res += fk[i]*(xk[i+1] - xk[i])
  res2h = 0
  for i in range(0, len(xk)-1, 2):
    res2h += fk[i]*(xk[i+2] - xk[i])
  return res + (res - res2h)

def trape_R(xk, fk):
  res = 0
  for i in range(len(xk)-1):
    res += (fk[i+1] + fk[i])/2*(xk[i+1] - xk[i])
  res2h=0
  for i in range(0, len(xk)-1, 2):
    res2h += (fk[i+2] + fk[i])/2*(xk[i+2] - xk[i])
  return res + (res - res2h)/3

#================ Метод прямоугольников =========================

I = 0
eps = list()
grid = list()
for i in range(1, 12):
  x = np.linspace(0, 5, pow(2, i))
  f = func(x)
  I = rect(x, f)
  eps.append(np.abs(I - ans))
  grid.append(x[1])
print('rect = ', I)

plt.figure(1, figsize=(5, 5))
plt.autoscale(tight=True)
plt.ylabel(r'$\ln(\delta)$')
plt.xlabel(r'$\ln(h)$')
plt.title("Метод прямоугольников")
E_x = np.log(grid)
E_y = np.log(eps)
plt.plot(E_x, E_y, linewidth=1, color='blue')
plt.errorbar(E_x, E_y, yerr=0.00, xerr=0.00, fmt='.', ecolor='red', color='red')
print("Порядок сходимости: ", (E_y[-1] - E_y[2])/(E_x[-1] - E_x[2]))
plt.grid()

#================ Метод трапеций =========================

print()
I = 0
eps = list()
grid = list()
for i in range(1, 12):
  x = np.linspace(0, 5, pow(2, i))
  f = func(x)
  I = trape(x, f)
  eps.append(np.abs(I - ans))
  grid.append(x[1])
print('trape = ', I)

plt.figure(2, figsize=(5, 5))
plt.autoscale(tight=True)
plt.ylabel(r'$\ln(\delta)$')
plt.xlabel(r'$\ln(h)$')
plt.title("Метод трапеций")
E_x = np.log(grid)
E_y = np.log(eps)
plt.plot(E_x, E_y, linewidth=1, color='blue')
plt.errorbar(E_x, E_y, yerr=0.00, xerr=0.00, fmt='.', ecolor='red', color='red')
print("Порядок сходимости: ", (E_y[-1] - E_y[2])/(E_x[-1] - E_x[2]))
plt.grid()

#================ Метод Симпсона =========================

print()
I = 0
eps = list()
grid = list()
for i in range(1, 12):
  x = np.linspace(0, 5, pow(2, i)+1)
  f = func(x)
  I = simp(x, f)
  eps.append(np.abs(I - ans))
  grid.append(2*x[1])
print('simp = ', I)

plt.figure(3, figsize=(5, 5))
plt.autoscale(tight=True)
plt.ylabel(r'$\ln(\delta)$')
plt.xlabel(r'$\ln(h)$')
plt.title("Метод Симпсона")
E_x = np.log(grid)
E_y = np.log(eps)
plt.plot(E_x, E_y, linewidth=1, color='blue')
plt.errorbar(E_x, E_y, yerr=0.00, xerr=0.00, fmt='.', ecolor='red', color='red')
print("Порядок сходимости: ", (E_y[-1] - E_y[3])/(E_x[-1] - E_x[3]))
plt.grid()

#================ Метод Рунге прямоугольников =========================

print()
I = 0
eps = list()
grid = list()
for i in range(1, 12):
  x = np.linspace(0, 5, pow(2, i)+1)
  f = func(x)
  I = rect_R(x, f)
  eps.append(np.abs(I - ans))
  grid.append(x[1])
print('rect_R = ', I)

plt.figure(4, figsize=(5, 5))
plt.autoscale(tight=True)
plt.ylabel(r'$\ln(\delta)$')
plt.xlabel(r'$\ln(h)$')
plt.title("Метод прямоугольников + Поправка")
E_x = np.log(grid)
E_y = np.log(eps)
plt.plot(E_x, E_y, linewidth=1, color='blue')
plt.errorbar(E_x, E_y, yerr=0.00, xerr=0.00, fmt='.', ecolor='red', color='red')
print("Порядок сходимости: ", (E_y[-1] - E_y[2])/(E_x[-1] - E_x[2]))
plt.grid()

#================ Метод трапеций Рунге =========================

print()
I = 0
eps = list()
grid = list()
for i in range(1, 12):
  x = np.linspace(0, 5, pow(2, i)+1)
  f = func(x)
  I = trape_R(x, f)
  eps.append(np.abs(I - ans))
  grid.append(x[1])
print('trape_R = ', I)

plt.figure(5, figsize=(5, 5))
plt.autoscale(tight=True)
plt.ylabel(r'$\ln(\delta)$')
plt.xlabel(r'$\ln(h)$')
plt.title("Метод трапеций + поправка")
E_x = np.log(grid)
E_y = np.log(eps)
plt.plot(E_x, E_y, linewidth=1, color='blue')
plt.errorbar(E_x, E_y, yerr=0.00, xerr=0.00, fmt='.', ecolor='red', color='red')
print("Порядок сходимости: ", (E_y[-1] - E_y[2])/(E_x[-1] - E_x[2]))
plt.grid()

#================= Несобственный интеграл =======================

def func_nes(x):
  return np.exp(-pow(x, 2))

def monte_carlo(itNum):
  res = 0
  for i in range(itNum):
    gamma = rd.random()
    xi = -1/3*np.log(gamma)
    res += func_nes(xi)/(3*np.exp(-3*xi))
  res = res/itNum
  return res

print()
I = 0
eps = list()
grid = list()

I = monte_carlo(1000000)
eps.append(np.abs(I - np.sqrt(np.pi)/2))
grid.append(1/np.sqrt(i))
print('Monte Carlo = ', I)

print("Погрешность: ", eps[-1])





plt.show()