import numpy as np 
x = np.linspace(0, 3, 1000)

def func(arg):
  if isinstance(arg, list):
    return [np.sin(3*x)*np.exp(-x**2/2) for x in arg]
  return np.sin(3*arg)*np.exp(-arg**2/2)

def diff1(arg):
  if isinstance(arg, list):
    return [(3*np.cos(3*x)-x*np.sin(3*x))*np.exp(-x**2/2) for x in arg]
  return (3*np.cos(3*arg)-arg*np.sin(3*arg))*np.exp(-arg**2/2)

def diff2(arg):
  if isinstance(arg, list):
    return [((x**2-10)*np.sin(3*x) - 6*x*np.cos(3*x))*np.exp(-x**2/2) for x in arg]
  return ((arg**2-10)*np.sin(3*arg) - 6*x*np.cos(3*arg))*np.exp(-arg**2/2)

def diff1_p1(arg, func_grid):
  res = list()
  res.append(0)
  h = arg[1]
  for i in range(1, len(arg)):
    res.append((func_grid[i] - func_grid[i-1])/h)
  res[0] = res[1]
  return res

def diff1_p2(arg, func_grid):
  res = list()
  res.append(0)
  h=arg[1]
  for i in range(1, len(arg) - 1):
    res.append((func_grid[i+1] - func_grid[i-1])/(2*h))
  res[0] = (-3*func_grid[0] + 4*func_grid[1] - func_grid[2])/(2*h)
  res.append((3*func_grid[-1] - 4*func_grid[-2] + func_grid[-3])/(2*h))
  return res

def diff2_p1(arg, func_grid):
  res = list()
  res.append(0)
  h=arg[1]
  for i in range(1, len(arg) - 1):
    res.append((func_grid[i-1] - 2*func_grid[i] + func_grid[i+1])/(h**2))
  res[0] = res[1]
  res.append(res[-1])
  return res

def diff2_p2(arg, func_grid):
  res = list()
  res.append(0)
  h=arg[1]
  for i in range(1, len(arg) - 1):
    res.append((func_grid[i-1] - 2*func_grid[i] + func_grid[i+1])/(h**2))
  res[0] = (2*func_grid[0] - 5*func_grid[1] + 4*func_grid[2] - func_grid[3])/(h**2)
  res.append((2*func_grid[-1] - 5*func_grid[-2] + 4*func_grid[-3] - func_grid[-4])/(h**2))
  return res

def countE(fu_lst, gr_fu_lst):
  E = [np.abs(fu_lst[i] - gr_fu_lst[i]) for i in range(len(fu_lst))]
  return max(E)