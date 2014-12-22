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
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

# Model parameters

a = (1./np.sqrt(3)) * 1.3 # Condition for the fix point to be in the outer branch a/e > 1/ sqrt(3)
e = 0.05 # Condition for excitable system e << 1

dt = 0.1

start  = 0.0
finish = 10000.0

def main():
		
	# Initial conditions:
	u_init = -a * 0.4
	v_init = -a + a**3
	
	plt.figure()
	ax = plt.gca()
	plt.ylim([-2,2])
	plt.xlim([0,50])
	plt.xlabel(r'$t$')
	plt.ylabel(r'$u$')
	

	t_plot, u_plot = evolve_FHN(u_init, v_init)
	plot_FHN(t_plot, u_plot, 'b', u_init)
	u_init = -a * 0.8
	t_plot, u_plot = evolve_FHN(u_init, v_init)
	plot_FHN(t_plot, u_plot, 'm', u_init)
	plt.axhline(y=-a,linewidth=2, color = 'g', linestyle="--")
	plt.legend(loc='best')
	plt.text(1.,1.,r'u_{0} = '+str(-a), transform=ax.transAxes, horizontalalignment='right', verticalalignment='bottom')
	plt.show()

def plot_FHN(t_plot, u_plot, line_color, u_init):
	plt.plot(t_plot, u_plot, line_color , linewidth=3,label=r'$u(t=0) = $'+str(u_init))
	
	

def step_u(u_old, v_old, dt):
	"""
	Implement the step of equation \dot(u) = u - u^3 - v
	"""
	return u_old + dt*(u_old - u_old**3 - v_old)

def step_v(u_old, v_old, dt):
	"""
	Implement the step of equation \dot(v) = e(u + a)
	"""
	return v_old + dt*(e*(u_old + a))
	

def evolve_FHN(u_init, v_init):
	
	t_plot = []
	u_plot = []
	
	t_plot.append(start)
	u_plot.append(u_init)
	
	u_old = u_init
	v_old = v_init	
	
		
	for t in np.arange(start+dt, finish+dt, dt):
		u_new = step_u(u_old, v_old, dt)
		v_new = step_v(u_old, v_old, dt)
		u_old = u_new
		v_old = v_new
		t_plot.append(t)
		u_plot.append(u_new)
		
	return t_plot, u_plot


if __name__ == "__main__":
	main()
