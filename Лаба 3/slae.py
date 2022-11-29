import numpy as np 
import random as rd

def Euc(a, b):
  dim = a.shape[0]
  s = 0
  for i in range(dim):
    s += a[i]*b[i]
  return s

def Hauss(B, fb, dim):
  for j in range(dim):
    maxInd = j
    for i in range(j, dim):
      if np.abs(B[i, j]) > np.abs(B[maxInd, j]):
        maxInd = i
    buff = np.array(B[j, :])
    B[j, :] = B[maxInd, :]
    B[maxInd, :] = buff
    buff = fb[j]
    fb[j] = fb[maxInd]
    fb[maxInd] = buff
    denom = B[j, j]
    B[j, :] /= denom
    fb[j] /= denom
    if(j < dim):
      for i in range(j+1, dim):
        a = B[i,j]
        fb[i] -= fb[j]*a
        B[i, :] -= B[j, :]*a
        for y in range(j):
          B[i, y] = 0
  xb = np.zeros((dim, 1))
  for i in range(dim-1, -1, -1):
    xb[i] = (fb[i] - B[i, :].dot(xb))
  return xb

def Through(B, f2, dim):
  koef = np.zeros((2, dim))
  koef[0, 0] = -B[0, 1]/B[0, 0]
  koef[1, 0] = f2[0]/B[0, 0]
  for i in range(1, dim):
    y = B[i, i] + B[i, i-1]*koef[0, i-1]
    if i < dim-1:
      koef[0, i] = -B[i, i+1]/y
    koef[1, i] = (f2[i] - koef[1, i-1]*B[i, i-1])/y
  xb = np.zeros((dim,1), dtype=np.float64)
  xb[dim-1] = koef[1, dim-1]
  for i in range(dim-2, -1, -1):
    xb[i] = koef[0,i]*xb[i+1] + koef[1, i]
  return xb

#=============== Метод Якоби ================

def Yakobi(B, fb, dim, it):
  xb = np.zeros((dim,1), dtype=np.float64)
  for i in range(dim):
    xb[i, 0] = rd.random()*1000
  for n in range(it):
    xn = np.array(xb)
    for i in range(dim): 
      xb[i] = -np.sum(B[i, :].dot(xn)/B[i, i]) + xn[i] + fb[i]/B[i, i]
  return(xb)

#================== Метод Зейделя =================

def Zeidel(B, fb, dim, it):
  xb = np.zeros((dim,1), dtype=np.float64)
  for i in range(dim):
    xb[i, 0] = rd.random()*1000
  for n in range(it):
    for i in range(dim): 
      xb[i] = 0
      xb[i] = -np.sum(B[i, :].dot(xb)/B[i, i]) + fb[i]/B[i, i]
  return xb

#============= Метод наискорейшего спуска ================

def FstDesc(B, fb, dim, itNum):
  xn = np.zeros((dim,1), dtype=np.float64)
  for i in range(dim):
    xn[i, 0] = rd.random()*1000
  for n in range(itNum):
    rn = B.dot(xn) - fb
    tn = Euc(rn, rn)/Euc(B.dot(rn), rn)
    xn = xn - tn*rn
  return xn

#============ Метод минимальной невязки ==================

def MinRes(B, fb, dim, itNum):
  xn = np.zeros((dim,1), dtype=np.float64)
  for i in range(dim):
    xn[i, 0] = rd.random()
  for n in range(itNum):
    rn = B.dot(xn) - fb
    tn = Euc(rn, B.dot(rn))/Euc(B.dot(rn), B.dot(rn))
    xn = xn - tn*rn
  return xn