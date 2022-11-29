import numpy as np 
import matplotlib.pyplot as plt

seg = np.array([0, 3.0])
def func1(x):
  return np.cos(x)
def diff1(x):
  return -np.sin(x)
def func2(x):
  return x**2
def diff2(x):
  return 2*x

def Func(x, y):
  return np.array([[x*x + y*y - 5], [x + y - 3]], dtype=np.float64)
def Fdiff(x, y):
  return np.array([[1, -2*y], [-1, 2*x]], dtype=np.float64)/(2*(x-y))

#------------Half division-------------
err = [0.1/2**i for i in range(50)]
iterations = list()
for eps in err:
  iteration = 1
  seg = np.array([0, 3.0])
  while (True):
    mid = (seg[0] + seg[1])/2
    if func1(mid) > eps:
      seg[0] = mid
    elif func1(mid) < -eps:
      seg[1] = mid
    else:
      break
    iteration += 1
  iterations.append(iteration)

plt.figure(1, figsize=(6, 6))
plt.autoscale(tight=True)
plt.xlabel(r'$\log(\varepsilon)$')
plt.ylabel('iterations')
plt.title(r'$f(x) = \cos(x)$')
plt.plot(np.log(err), iterations, linewidth=2)
plt.grid()

  #-----------Newton for eq----------
delta = list()
for n in range(50):
  x = np.array([2.0, 2, 2])
  for i in range(n):
    x[0] = x[1]
    x[1] = x[2]
    x[2] = x[1] - func2(x[1])/diff2(x[1])
  delta.append(x[2])
q = x[2]/x[1]
print("q = ",q)

delta2 = list()
for n in range(50):
  x = np.array([2.0, 2, 2])
  koef = 1/diff2(x[0])
  for i in range(n):
    x[0] = x[1]
    x[1] = x[2]
    x[2] = x[1] - func2(x[1])*koef
  delta2.append(x[2])

plt.figure(2, figsize=(12, 6))
plt.subplot(121)
plt.autoscale(tight=True)
plt.ylabel(r'$\Delta$')
plt.xlabel('iterations')
plt.title(r'$f(x) = x^2, Newton$')
plt.plot(range(50), delta, linewidth=2)
plt.grid()
plt.subplot(122)
plt.autoscale(tight=True)
plt.ylabel(r'$\Delta$')
plt.xlabel('iterations')
plt.title(r'$f(x) = x^2$')
plt.plot(range(50), delta2, linewidth=2)
plt.grid()

#------------Newton for sys------------
delta = list()
for n in range(50):
  x = np.array([[2], [3]], dtype=np.float64)
  for i in range(n):
    x = x - Fdiff(x[0,0], x[1,0]).dot(Func(x[0,0], x[1,0]))
  delta.append(np.sqrt((x[0,0] - 1)**2 + (x[1,0] - 2)**2))
print("x = ", x)
print("Delta_50 = ", delta[-1])
plt.figure(3, figsize=(6, 6))
plt.autoscale(tight=True)
plt.ylabel(r'$\Delta$')
plt.xlabel('iterations')
plt.title(r'$x^2+y^2=5, x+y =3$')
plt.plot(range(50), delta, linewidth=2)
plt.grid()


plt.show()
