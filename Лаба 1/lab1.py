import numpy as np 
import matplotlib.pyplot as plt

from diffLib import *

k_max = 17

#============ Function ===================

plt.figure(1, figsize=(6, 6))
plt.autoscale(tight=True)
plt.xlabel(r'$x$')
plt.ylabel(r'$f(x)$')
plt.title(r'$f(x) = \sin(3x)*\exp(-x^2/2)$')
plt.plot(x, func(x), linewidth=2)
plt.grid()

#=========== First diff ==================
#first pow

plt.figure(2, figsize=(10, 10))
plt.subplot(221)
plt.title(r"$f', n = 1$")
plt.autoscale(tight=True)
plt.plot(x, diff1(x), linewidth=2)
plt.grid()

plt.subplot(222)
k=2
x_grid = [i/2**k for i in range((3)*2**k+1)]  #making grid with h = 1/2
plt.plot(x, diff1(x), linewidth=2)
plt.plot(x_grid, diff1_p1(x_grid, func(x_grid)), linewidth=1, color='red')
plt.grid()
plt.title(r'$h = 2^{-2}$')

plt.subplot(223)
k=6
x_grid = [i/2**k for i in range((3)*2**k+1)]  #making grid with h = 1/64
plt.plot(x, diff1(x), linewidth=2)
plt.plot(x_grid, diff1_p1(x_grid, func(x_grid)), linewidth=1, color='red')
plt.grid()
plt.title(r'$h = 2^{-6}$')

plt.subplot(224)
Epsilon = [list(), list()]
for k in range(1, k_max):
  print(k)
  x_grid = [i/2**k for i in range(3*2**k+1)]  #making grid with h = 1/2^k
  Epsilon[0].append(x_grid[1])
  Epsilon[1].append(countE(diff1(x_grid), diff1_p1(x_grid, func(x_grid))))
E_x = [np.log(x) for x in Epsilon[0]]
E_y = [np.log(x) for x in Epsilon[1]]
plt.errorbar(E_x, E_y, yerr=0.00, xerr=0.00, fmt='.', ecolor='red', color='red')
plt.plot(E_x, E_y, linewidth=1, color='red')
print((E_y[-1] - E_y[1])/(E_x[-1] - E_x[1]))
plt.grid()
plt.title(r'$tg = 0,997$')

#second pow

plt.figure(3, figsize=(10, 10))
plt.subplot(221)
plt.title(r"$f', n = 2$")
plt.autoscale(tight=True)
plt.plot(x, diff1(x), linewidth=2)
plt.grid()

plt.subplot(222)
k=2
x_grid = [i/2**k for i in range((3)*2**k+1)]  #making grid with h = 1/2
plt.plot(x, diff1(x), linewidth=2)
plt.plot(x_grid, diff1_p2(x_grid, func(x_grid)), linewidth=1, color='red')
plt.grid()
plt.title(r'$h = 2^{-2}$')

plt.subplot(223)
k=6
x_grid = [i/2**k for i in range((3)*2**k+1)]  #making grid with h = 1/64
plt.plot(x, diff1(x), linewidth=2)
plt.plot(x_grid, diff1_p2(x_grid, func(x_grid)), linewidth=1, color='red')
plt.grid()
plt.title(r'$h = 2^{-6}$')

plt.subplot(224)
Epsilon = [list(), list()]
for k in range(1, k_max):
  print(k)
  x_grid = [i/2**k for i in range(3*2**k+1)]  #making grid with h = 1/2^k
  Epsilon[0].append(x_grid[1])
  Epsilon[1].append(countE(diff1(x_grid), diff1_p2(x_grid, func(x_grid))))
E_x = [np.log(x) for x in Epsilon[0]]
E_y = [np.log(x) for x in Epsilon[1]]
plt.errorbar(E_x, E_y, yerr=0.00, xerr=0.00, fmt='.', ecolor='red', color='red')
plt.plot(E_x, E_y, linewidth=1, color='red')
print((E_y[-1] - E_y[2])/(E_x[-1] - E_x[2]))
plt.grid()
plt.title(r'$tg = 1,988$')

#=========== Second diff ==================

#first pow

plt.figure(4, figsize=(10, 10))
plt.subplot(221)
plt.title(r"$f'', n = 1$")
plt.autoscale(tight=True)
plt.plot(x, diff2(x), linewidth=2)
plt.grid()

plt.subplot(222)
k=2
x_grid = [i/2**k for i in range((3)*2**k+1)] 
plt.plot(x, diff2(x), linewidth=2)
plt.plot(x_grid, diff2_p1(x_grid, func(x_grid)), linewidth=1, color='red')
plt.grid()
plt.title(r'$h = 2^{-2}$')

plt.subplot(223)
k=6
x_grid = [i/2**k for i in range((3)*2**k+1)]
plt.plot(x, diff2(x), linewidth=2)
plt.plot(x_grid, diff2_p1(x_grid, func(x_grid)), linewidth=1, color='red')
plt.grid()
plt.title(r'$h = 2^{-6}$')

plt.subplot(224)
Epsilon = [list(), list()]
for k in range(1, k_max):
  print(k)
  x_grid = [i/2**k for i in range(3*2**k+1)]
  Epsilon[0].append(x_grid[1])
  Epsilon[1].append(countE(diff2(x_grid), diff2_p1(x_grid, func(x_grid))))
E_x = [np.log(x) for x in Epsilon[0]]
E_y = [np.log(x) for x in Epsilon[1]]
plt.errorbar(E_x, E_y, yerr=0.00, xerr=0.00, fmt='.', ecolor='red', color='red')
plt.plot(E_x, E_y, linewidth=1, color='red')
print((E_y[-1] - E_y[2])/(E_x[-1] - E_x[2]))
plt.grid()
plt.title(r'$tg = 0,992$')

#second pow

plt.figure(5, figsize=(10, 10))
plt.subplot(221)
plt.title(r"$f'', n = 2$")
plt.autoscale(tight=True)
plt.plot(x, diff2(x), linewidth=2)
plt.grid()

plt.subplot(222)
k=2
x_grid = [i/2**k for i in range((3)*2**k+1)]
plt.plot(x, diff2(x), linewidth=2)
plt.plot(x_grid, diff2_p2(x_grid, func(x_grid)), linewidth=1, color='red')
plt.grid()
plt.title(r'$h = 2^{-2}$')

plt.subplot(223)
k=6
x_grid = [i/2**k for i in range((3)*2**k+1)]
plt.plot(x, diff2(x), linewidth=2)
plt.plot(x_grid, diff2_p2(x_grid, func(x_grid)), linewidth=1, color='red')
plt.grid()
plt.title(r'$h = 2^{-6}$')

plt.subplot(224)
Epsilon = [list(), list()]
for k in range(1, k_max):
  print(k)
  x_grid = [i/2**k for i in range(3*2**k+1)]
  Epsilon[0].append(x_grid[1])
  Epsilon[1].append(countE(diff2(x_grid), diff2_p2(x_grid, func(x_grid))))
E_x = [np.log(x) for x in Epsilon[0]]
E_y = [np.log(x) for x in Epsilon[1]]
plt.errorbar(E_x, E_y, yerr=0.00, xerr=0.00, fmt='.', ecolor='red', color='red')
plt.plot(E_x, E_y, linewidth=1, color='red')
koef = (E_y[12] - E_y[5])/(E_x[12] - E_x[5])
print(koef)
xlin = np.linspace(E_x[5], E_x[12], 1000)
plt.plot(xlin, koef*xlin + E_x[5] + 6.3)
plt.grid()
plt.title(r'$tg = 1,966$')


plt.show()