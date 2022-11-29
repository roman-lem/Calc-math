import numpy as np 
import matplotlib.pyplot as plt
import random as rd
from slae import *

dim = 100

def rhoEuc(x, xb, dim):
  s = 0
  for i in range(dim):
    s += (x[i] - xb[i])**2
  return np.sqrt(s)


A = np.zeros((dim, dim), dtype=np.float64)
x = np.zeros((dim,1), dtype=np.float64)
for i in range(dim):
  x[i, 0] = 100


rd.seed(version = 2)

for i in range(dim):
  for j in range(dim):
    if i != j:
      A[i, j] = rd.uniform(-100000, 100000)
      A[i, i] += np.abs(A[i][j]) + rd.random()*1000

f = A.dot(x)

#print(A)
#print(x)
#print(f)
#print("="*40)

#===================== Метод гаусса =========================

B = np.array(A)
fb = np.array(f)
xb = Hauss(B, fb, dim)

print("Метод Гаусса")
print("rhoEuc(x, x_*) = ", rhoEuc(x, xb, dim))

#=============== Метод прогонки =============

B = np.zeros((dim, dim))
B[0, 0] = rd.uniform(-100000, 100000)
B[0, 1] = rd.uniform(-100000, 100000)
for i in range(1, dim-1):
  B[i, i-1] = rd.uniform(-100000, 100000)
  B[i,  i ] = rd.uniform(-100000, 100000)
  B[i, i+1] = rd.uniform(-100000, 100000)
B[dim-1, dim-2] = rd.uniform(-100000, 100000)
B[dim-1, dim-1] = rd.uniform(-100000, 100000)

f2 = B.dot(x)

xb = Through(B, f2, dim)
print()
print("Метод прогонки")
print("rhoEuc(x, x_*) = ", rhoEuc(x, xb, dim))


#=============== Метод Якоби ================

B = np.array(A)
fb = np.array(f)

delta = list()
for itNum in range(20):
  xb = Yakobi(B, fb, dim, itNum)
  delta.append(rhoEuc(x, xb, dim))


plt.figure(1, figsize=(12, 6))
plt.subplot(121)
plt.autoscale(tight=True)
plt.ylabel(r'$\Delta$')
plt.xlabel('iterations')
plt.title('Метод Якоби')
plt.plot(range(1,21), delta, linewidth=2)
plt.grid()
print()
print("Метод Якоби")
print('Delta(5) = ', delta[4])
print('Delta(20) = ', delta[19])

#================== Метод Зейделя =================

B = np.array(A)
fb = np.array(f)

delta = list()
for itNum in range(20):
  xb = Zeidel(B, fb, dim, itNum)
  delta.append(rhoEuc(x, xb, dim))


plt.subplot(122)
plt.autoscale(tight=True)
plt.ylabel(r'$\Delta$')
plt.xlabel('iterations')
plt.title('Метод Зейделя')
plt.plot(range(1,21), delta, linewidth=2)
plt.grid()
print()
print("Метод Зейделя")
print('Delta(5) = ', delta[4])
print('Delta(20) = ', delta[19])

#============= Метод наискорейшего спуска ================

B = np.array(A)
B = B.transpose().dot(B)
fb = A.transpose().dot(f)


delta = list()
for itNum in range(20):
  xn = FstDesc(B, fb, dim, itNum)
  delta.append(rhoEuc(x, xn, dim))

plt.figure(2, figsize=(12, 6))
plt.subplot(121)
plt.autoscale(tight=True)
plt.ylabel(r'$\Delta$')
plt.xlabel('iterations')
plt.title('Метод наискорейшего спуска')
plt.plot(range(1,21), delta, linewidth=2)
plt.grid()
print()
print("Метод наискорейшего спуска")
print('Delta(5) = ', delta[4])
print('Delta(20) = ', delta[19])

#============ Метод минимальной невязки ==================

B = np.array(A)
B = B.transpose().dot(B)
fb = A.transpose().dot(f)

delta = list()
for itNum in range(20):
  xn = MinRes(B, fb, dim, itNum)
  delta.append(rhoEuc(x, xn, dim))

plt.figure(2, figsize=(12, 6))
plt.subplot(122)
plt.autoscale(tight=True)
plt.ylabel(r'$\Delta$')
plt.xlabel('iterations')
plt.title('Метод минимальной невязки')
plt.plot(range(1,21), delta, linewidth=2)
plt.grid()
print()
print("Метод минимальной невязки")
print('Delta(5) = ', delta[4])
print('Delta(20) = ', delta[19])



plt.show()
