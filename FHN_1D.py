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

a0 = (1./np.sqrt(3)) * 1.3 # Condition for the fix point to be in the outer branch a/e > 1/ sqrt(3)
a1 = 3.
e = 0.02 # Condition for excitable system e << 1
delta = 0.1 # Diffusion coefficient

start  = 0.0
finish = .0005

L = 10.
N = 100.
l = np.linspace(0,L,N)

dx = L/N
dx2 = dx**2
dt = dx2 / 100


	
fig, ax = plt.subplots()

	
# Initial conditions:
u_init = 0.2*((np.random.rand(N))-0.5) - a0   # random initial conditions around u_0
v_init = 0.2*((np.random.rand(N))-0.5) -a0 + a0**3 # random initial conditions around v_0

line, = ax.plot(l, u_init)

u = u_init.copy()
v = v_init.copy()

def updatefig(*args):
    global l,u,v,dt,dx2
    u += dt*(u - u**3 - v + laplacian(u, dx2))
    v += dt*dt*(e*(u - a1*v + a0) + delta*laplacian(v,dx2))
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
	return (np.roll(var,1,axis=0) + np.roll(var,-1,axis=0) -2*var)/dh2


def step_FHN(u_old, v_old, dt, dx2):
	"""
	
	"""
	return step_u(u_old, v_old, dt, dx2), step_v(u_old, v_old, dt, dx2)
	
def evolve_FHN(u_init, v_init):
		
	
	t_plot = []
	u_plot = []
	
	t_plot.append(start)
	u_plot.append(u_init)
	
	u_old = u_init
	v_old = v_init	
	
		
	for t in np.arange(start+dt, finish+dt, dt):
		u_new = step_u(u_old, v_old, dt, dx2)
		v_new = step_v(u_old, v_old, dt, dx2)
		u_old = u_new
		v_old = v_new
		t_plot.append(t)
		u_plot.append(u_new)
		
	return t_plot, u_plot


ani = animation.FuncAnimation(fig, updatefig, interval=100, blit=False)
plt.show()
