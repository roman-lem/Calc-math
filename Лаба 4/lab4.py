import numpy as np 
import matplotlib.pyplot as plt
from slae import *

class My_Poly():

  def get_coef(self):
    return np.array(self.coef)

  def count(self, x):
    ret=0
    for i in range(len(self.coef)):
      ret += self.coef[i]*pow(x, i)
    return ret

  def mulxn(self, n):
    res = np.insert(self.coef, 0, np.zeros(n))
    return self.__class__(res)

  def mulint(self, a):
    res = self.coef*a
    return self.__class__(res)

  def getPow():
    return len(self.coef)-1

  def __init__(self, coef):
    self.coef = np.array(coef)

  def __call__(self, arg):
    if isinstance(arg, (int, float)):
      return self.count(arg)
    elif isinstance(arg, (np.ndarray, list)):
      ret = np.zeros(len(arg))
      for i in range(len(arg)):
        ret[i] = self.count(arg[i])
      return ret

  def __add__(self, other):
    othercoef = other.get_coef()
    lendiff = len(othercoef) - len(self.get_coef())
    if (lendiff > 0):
      retcoef = othercoef + np.append(self.get_coef(), np.zeros(lendiff))
    elif (lendiff < 0):
      retcoef = np.append(othercoef, np.zeros(-lendiff)) + self.get_coef()
    elif (lendiff == 0):
      retcoef = othercoef + self.get_coef()
    return self.__class__(retcoef)

  def __mul__(self, other):
    if isinstance(other, (int, float)):
      return self.mulint(other)
    elif isinstance(other, self.__class__):
      othercoef = other.get_coef()
      ret = self.mulint(othercoef[0])
      for i in range(1, len(othercoef)):
        ret = ret + self.mulint(othercoef[i]).mulxn(i)
      return ret

  def __str__(self):
    coefs = self.get_coef()
    s = str(coefs[0])
    for i in range(1, len(coefs)):
      if coefs[i] < 0:
        s += ' - ' + str(-coefs[i])+ 'x^'+str(i)
      if coefs[i] > 0:
        s += ' + ' + str(coefs[i])+ 'x^'+str(i)
    return s

  def __truediv__(self, other):
    if other != 0:
      coefs = self.get_coef()/other
      return self.__class__(coefs)

class Lagrange(My_Poly):

  def __init__(self, xk, fk):
    self.coef = np.array([0])
    for i in range(len(xk)):
      l = My_Poly([1])
      for j in range(len(xk)):
        if j != i:
          l = l*My_Poly([-xk[j], 1])/(xk[i] - xk[j])
      self.coef =  (My_Poly(self.get_coef()) + l*fk[i]).get_coef()

class Newton(My_Poly):

  def divDiff(self, xk, fk, n, ini):
    if n == 1:
      return (fk[ini+n] - fk[ini])/(xk[ini+n] - xk[ini])
    else :
      return (self.divDiff(xk, fk, n-1, ini+1) - self.divDiff(xk, fk, n-1, ini))/(xk[ini+n] - xk[ini])

  def __init__(self, xk, fk):
    self.x = xk
    self.f = fk
    self.coef = np.array([fk[0]])
    poly = My_Poly([1])
    for i in range(len(xk)-1):
      poly *= My_Poly([-xk[i], 1])
      self.coef =  (My_Poly(self.get_coef()) + poly*self.divDiff(xk, fk, i+1, 0)).get_coef()

  def append(self, x, f):
    self.x.append(x)
    self.f.append(f)
    poly = My_Poly([1])
    for i in range(self.getPow()-1):
      poly *= My_Poly([-self.x[i], 1])
    self.coef =  (My_Poly(self.get_coef()) + poly*self.divDiff(xk, fk, self.getPow(), 0)).get_coef()

class CubeSpline():
  def get_parts(self):
    return np.array(self.parts)

  def count(self, x):
    if x < self.x[0] or x > self.x[-1]:
      return np.NaN
    ind = 1 
    while(self.x[ind] < x):
      ind += 1
    ind -= 1
    poly = self.parts[ind]
    return poly(x)

  def __init__(self, xk, fk, dfSide):
    self.x = np.array(xk)
    N = len(xk)-1
    A = np.zeros((4*N, 4*N))
    fb = np.zeros((4*N, 1))
    for i in range(4):
      A[0][i] = pow(xk[0], 3-i)
      fb[0] = fk[0]
    for i in range(3):
      A[1][i] = (3 - i)*pow(xk[0], 2-i)
      fb[1] = dfSide[0]
    for i in range(4):
      A[2][(N-1)*4+i] = pow(xk[N], 3-i)
      fb[2] = fk[N]
    for i in range(3):
      A[3][(N-1)*4+i] = (3 - i)*pow(xk[N], 2-i)
      fb[3] = dfSide[1]
    for i in range(1, N):
      for j in range(4):
        A[i*4][i*4+j] = pow(xk[i], 3-j)
        fb[i*4] = fk[i]
      for j in range(4):
        A[i*4+1][(i-1)*4+j] = pow(xk[i], 3-j)
        fb[i*4+1] = fk[i]
      for j in range(3):
        A[i*4+2][(i-1)*4+j] = -(3-j)*pow(xk[i], 2-j)
      for j in range(3):
        A[i*4+2][i*4+j] = (3-j)*pow(xk[i], 2-j)
      fb[i*4+2] = 0
      A[i*4+3][(i-1)*4] = -6*xk[i]
      A[i*4+3][(i-1)*4+1] = -2*xk[i]
      A[i*4+3][i*4] = 6*xk[i]
      A[i*4+3][i*4+1] = 2*xk[i]
      fb[i*4+3] = 0
    forPoly = Hauss(np.array(A), np.array(fb), 4*N)
    B = A.transpose().dot(A)
    fb = A.transpose().dot(fb)
    forPoly = MinRes(B, fb, 4*N, 100, forPoly)
    self.parts = list()
    for i in range(N):
      self.parts.append(My_Poly([forPoly[i*4 + j] for j in range(3, -1, -1)]))

  def __call__(self, arg):
    if isinstance(arg, (int, float)):
      return self.count(arg)
    elif isinstance(arg, (np.ndarray, list)):
      ret = np.zeros(len(arg))
      for i in range(len(arg)):
        ret[i] = self.count(arg[i])
      return ret

  def __str__(self):
    N = len(xk)-1
    s = str(self.parts[0])
    for i in range(N):
      s += "Poly " + str(i) + ": " + str(self.parts[i])+'\n'
    return s



class My_Func():

  def __call__(self, arg):
    if isinstance(arg, (list, np.ndarray)):
      return np.array([np.sin(3*x)*np.exp(-x**2/2) for x in arg])
    return np.sin(3*arg)*np.exp(-arg**2/2)

class My_Diff():
  def __call__(self, arg):
    if isinstance(arg, (list, np.ndarray)):
      return np.array([(3*np.cos(3*x)-x*np.sin(3*x))*np.exp(-x**2/2) for x in arg])
    else:
      return (3*np.cos(3*arg)-arg*np.sin(3*arg))*np.exp(-arg**2/2)

#============ Интерполяция полиномом Лагранжа =================

f = My_Func()

x=np.linspace(0, 3, 10000)
xk = np.linspace(0, 3, 1001)
fk = f(xk)

plt.figure(1, figsize=(10, 10))
plt.subplot(221)
L = Lagrange(xk[::200], fk[::200])
plt.autoscale(tight=True)
plt.ylabel(r'$f(x)$')
plt.title(r'$f(x) = \sin(3x)*\exp(-x^2/2)$, 6т.')
plt.plot(x, f(x), linewidth=2)
plt.plot(x, L(x), linewidth=2)
plt.grid()

plt.subplot(222)
L = Lagrange(xk[::100], fk[::100])
plt.autoscale(tight=True)
plt.ylabel(r'$f(x)$')
plt.title(r'11т.')
plt.plot(x, f(x), linewidth=2)
plt.plot(x, L(x), linewidth=2)
plt.grid()

plt.subplot(223)
L = Lagrange(xk[::50], fk[::50])
plt.autoscale(tight=True)
plt.xlabel(r'$x$')
plt.ylabel(r'$f(x)$')
plt.title(r'21т.')
plt.plot(x, f(x), linewidth=2)
plt.plot(x, L(x), linewidth=2)
plt.grid()

epsilon = []
for i in range(2, 21):
  xk = np.linspace(0, 3, i)
  L = Lagrange(xk, f(xk))
  epsilon.append(max(np.abs(L(xk) - f(xk))))
plt.subplot(224)
plt.autoscale(tight=True)
plt.xlabel(r'$dot num$')
plt.ylabel(r'$\Delta$')
plt.title(r'optimal - 18т.')
plt.plot(range(2, 21), epsilon, linewidth=2)
plt.grid()

# #============ Интерполяция полиномом Ньютона =================

f = My_Func()

x = np.linspace(0, 3, 10000)
xk = np.linspace(0, 3, 1001)
fk = f(xk)

plt.figure(2, figsize=(10, 10))
plt.subplot(221)
N = Newton(xk[::200], fk[::200])
plt.autoscale(tight=True)
plt.ylabel(r'$f(x)$')
plt.title(r'$f(x) = \sin(3x)*\exp(-x^2/2)$, 6т.')
plt.plot(x, f(x), linewidth=2)
plt.plot(x, N(x), linewidth=2)
plt.grid()

plt.subplot(222)
N = Newton(xk[::100], fk[::100])
plt.autoscale(tight=True)
plt.ylabel(r'$f(x)$')
plt.title(r'11т.')
plt.plot(x, f(x), linewidth=2)
plt.plot(x, N(x), linewidth=2)
plt.grid()

plt.subplot(223)
N = Newton(xk[::50], fk[::50])
plt.autoscale(tight=True)
plt.xlabel(r'$x$')
plt.ylabel(r'$f(x)$')
plt.title(r'21т.')
plt.plot(x, f(x), linewidth=2)
plt.plot(x, N(x), linewidth=2)
plt.grid()

epsilon = []
for i in range(2, 25):
  xk = np.linspace(0, 3, i)
  N = Newton(xk, f(xk))
  epsilon.append(max(np.abs(N(xk) - f(xk))))
plt.subplot(224)
plt.autoscale(tight=True)
plt.xlabel(r'$dot num$')
plt.ylabel(r'$\Delta$')
plt.title(r'optimal - ?')
plt.plot(range(2, 25), epsilon, linewidth=2)
plt.grid()

#============ Интерполяция кубическим сплайном =================

f = My_Func()
df = My_Diff()

x = np.linspace(0, 3, 10000)
xk = np.linspace(0, 3, 1001)
fk = f(xk)

plt.figure(3, figsize=(10, 10))
plt.subplot(221)
xs = xk[::200]
fs = fk[::200]
Sp = CubeSpline(xs, fs, [df(xs[0]), df(xs[-1])])
plt.autoscale(tight=True)
plt.ylabel(r'$f(x)$')
plt.title(r'$f(x) = \sin(3x)*\exp(-x^2/2)$, 6т.')
plt.plot(x, f(x), linewidth=2)
plt.plot(x, Sp(x), linewidth=2)
plt.grid()

plt.subplot(222)
xs = xk[::100]
fs = fk[::100]
Sp = CubeSpline(xs, fs, [df(xs[0]), df(xs[-1])])
plt.autoscale(tight=True)
plt.ylabel(r'$f(x)$')
plt.title(r'11т.')
plt.plot(x, f(x), linewidth=2)
plt.plot(x, Sp(x), linewidth=2)
plt.grid()

plt.subplot(223)
xs = xk[::50]
fs = fk[::50]
Sp = CubeSpline(xs, fs, [df(xs[0]), df(xs[-1])])
plt.autoscale(tight=True)
plt.ylabel(r'$f(x)$')
plt.title(r'21т.')
plt.plot(x, f(x), linewidth=2)
plt.plot(x, Sp(x), linewidth=2)
plt.grid()

epsilon = []
for i in range(2, 25):
  xk = np.linspace(0, 3, i)
  Sp = CubeSpline(xk, f(xk), [df(xk[0]), df(xk[-1])])
  epsilon.append(max(np.abs(Sp(xk) - f(xk))))
plt.subplot(224)
plt.autoscale(tight=True)
plt.xlabel(r'$dot num$')
plt.ylabel(r'$\Delta$')
plt.title(r'optimal - ?')
plt.plot(range(2, 25), epsilon, linewidth=2)
plt.grid()



plt.show()
