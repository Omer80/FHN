#! /usr/bin/env python
"""
1D FHN integrator using explicit finite differences:
\dot(u) = u - u^3 - v
\dot(v) = e(u+a)
Periodic boundary conditions are assumed.
"""
__version__=1.0
__author__ = """Omer Tzuk (cliffon@gmail.com)"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

# Model parameters

a0 = 0.0 # Condition for the fix point to be in the outer branch a/e > 1/ sqrt(3)
a1 = 0.5
delta = 7.5 # Diffusion coefficient
e = delta / ((2.0-a1) + 2.0*np.sqrt(1.0-a1)) - 0.1 # Condition for excitable system e_c = delta / ((2.0-a) + 2.0*np.sqrt(1.0-a))

start  = 0.0
finish = .0005

L = 3.0
N = 1000
l = np.linspace(0,L,N)

dx = L/N
dx2 = dx**2
dt = dx2 / 100.0


# Random initial conditions:
u_init = 0.2*((np.random.rand(N))-0.5) - a0   # random initial conditions around u_0
v_init = 0.2*((np.random.rand(N))-0.5) -a0 + a0**3 # random initial conditions around v_0

# Sinus functions

#u_init = np.sin(np.linspace(0,L*2*np.pi,N))
#v_init = np.sin(np.linspace(0,L*2*np.pi,N))

u = u_init.copy()
v = v_init.copy()

def updatefig(*args):
    global l,u,v,dt,dx2
    u += dt*(u - u**3 - v + laplacian(u, dx2))
    v += dt*(e*(u - a1*v + a0) + delta*laplacian(v,dx2))
    line.set_ydata(u)
    return line,
	

	
def step_u(u_old, v_old, dt, dx2):
	"""
	Implement the step of equation \dot(u) = u - u^3 - v
	"""
	return u_old + dt*(u_old - u_old**3 - v_old + laplacian(u_old, dx2))

def step_v(u_old, v_old, dt, dx2):
	"""
	Implement the step of equation \dot(v) = e(u + a)
	"""
	return v_old + dt*(e*(u_old + a0) + delta*laplacian(v_old,dx2))
	
def laplacian(var, dh2):
	"""
	Implement the d^2(u)/dx^2 for periodic boundary condition
	"""
	numer = (3.0/2.0)*np.roll(var,1,axis=0) + (3.0/2.0)*np.roll(var,-1,axis=0) + (-49.0/18.0)*var
	numer += (-3.0/20.0)*np.roll(var,2,axis=0) + (-3.0/20.0)*np.roll(var,-2,axis=0)
	numer += (1.0/90.0)*np.roll(var,3,axis=0) + (1.0/90.0)*np.roll(var,-3,axis=0)
	return numer/dh2




if __name__=="__main__":
	fig, ax = plt.subplots()
	line, = ax.plot(l, u_init)
	ani = animation.FuncAnimation(fig, updatefig, interval=1, blit=False)
	plt.show()
