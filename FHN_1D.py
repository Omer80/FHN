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


par = {'a0' : 0.0, 
	   'a1' : 0.5,
	   'delta' : 7.5,
	   'e' : 2.0,
	   'L' : 3.0,
	   'N' : 1000,
	   }
	   
	   
def main():
	start  = 0.0
	step   = 1.0e-1
	finish = 100.0
	
	# Condition for excitable system e_c = delta / ((2.0-a) + 2.0*np.sqrt(1.0-a))
	ec = par['delta'] / ((2.0-par['a1']) + 2.0*np.sqrt(1.0-par['a1']))
	e =  ec - 0.1 
	l = np.linspace(0,par['L'],par['N'])
	dx = par['L']/par['N']
	dx2 = dx**2
	dt = 0.00001/(2.0 * dx**2)
	
	# Critical wavenumber 
	kc = np.sqrt((par['delta'] - ec*par['a1'])/ (2.0*par['delta']))
	
	par.update(dx=dx, dt=dt, e=e,dx2=dx2, ec=ec, kc=kc)
	
	
	# Random initial conditions:
	u_init = 0.2*((np.random.rand(par['N']))-0.5) - par['a0']   # random initial conditions around u_0
	v_init = 0.2*((np.random.rand(par['N']))-0.5) -par['a0'] + par['a0']**3 # random initial conditions around v_0
	
	# Sinus functions
	#u_init = np.sin(np.linspace(0,L*2*np.pi,N))
	#v_init = np.sin(np.linspace(0,L*2*np.pi,N))
	
	plt.figure()
	
	title=plt.title('time=%2.1f'%start)
	
	u_old = u_init.copy() 
	v_old = v_init.copy() 
	
	#plt.draw()
	plt.savefig('img_t_0.png', format='png', dpi=1000)
	image = 0
	t=start
	# start loop
	for tout in np.arange(start+step,finish+step,step):
	    while t < tout:
			#print "time:", t
			u_new = step_u(u_old, v_old, dt, dx2)
			v_new = step_v(u_old, v_old, dt, dx2)
			u_old = u_new
			v_old = v_new
			t+=dt
	    title.set_text('time=%2f'%(t))
	    plt.plot(l,u_new)
	    image +=1
	    plt.savefig('img_t_'+str(image)+'.png', format='png', dpi=1000)
	    #plt.draw()
	

	
def step_u(u_old, v_old, dt, dx2):
	"""
	Implement the step of equation \dot(u) = u - u^3 - v
	"""
	return u_old + dt*(u_old - u_old**3 - v_old + laplacian(u_old, dx2))

def step_v(u_old, v_old, dt, dx2):
	"""
	Implement the step of equation \dot(v) = e(u + a)
	"""
	return v_old + dt*(par['e']*(u_old + par['a0']) + par['delta']*laplacian(v_old,dx2))
	
def laplacian(var, dh2):
	"""
	Implement the d^2(u)/dx^2 for periodic boundary condition
	"""
	numer = (3.0/2.0)*np.roll(var,1,axis=0) + (3.0/2.0)*np.roll(var,-1,axis=0) + (-49.0/18.0)*var
	numer += (-3.0/20.0)*np.roll(var,2,axis=0) + (-3.0/20.0)*np.roll(var,-2,axis=0)
	numer += (1.0/90.0)*np.roll(var,3,axis=0) + (1.0/90.0)*np.roll(var,-3,axis=0)
	return numer/dh2



if __name__=="__main__":
	main()
