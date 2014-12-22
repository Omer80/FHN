#! /usr/bin/env python
"""
2D FHN integrator using explicit finite differences:
\dot(u) = u - u^3 - v + delta*laplacian
\dot(v) = e(u + a1*v +a0) + laplacian
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


par = {'a0': (1./np.sqrt(3)) * 0.5, 
       'a1':-2.0,
       'e':3.,
       'delta':50.,
       'n':(100,100),
       'l':(10,10),
       }

dx = float(par['l'][0]) / par['n'][0]
dy = float(par['l'][1]) / par['n'][1]
par.update(dx=dx, dy=dy)
X,Y=np.mgrid[0:par['n'][0],0:par['n'][1]]
X = X*dx
Y = Y*dy

dx2=dx**2 # To save CPU cycles, we'll compute Delta x^2
dy2=dy**2 # and Delta y^2 only once and store them.
h2 = dx2
h4 = h2**2

# initial condition
u0 = 0.5*(np.random.random((par['n'][0],par['n'][1]))- par['a0'])  # random initial conditions
v0 = 0.5*(np.random.random((par['n'][0],par['n'][1]))- par['a0'] + par['a0']**3) # random initial conditions
#u0 = 0.5* np.ones((par['n'][0],par['n'][1]))  # set initial conditions of constant value of 0.5

u0[1,:] = -par['a0'] + 1
start  = 0.0
step   = 0.1
finish = 1.0
#dt     = 0.1 # Semi-spectral time step
dt     = dx2 / 50000 # FDM time step


def main():
	
	          
	
	# plot first frame (t=start)
	plt.ion()
	plt.clf()
	ext = [0,par['l'][0],0,par['l'][1]]
	
	im=plt.imshow(u0.T,origin='lower', interpolation='nearest', extent=ext, cmap='jet')
	cbar=plt.colorbar()
	title=plt.title('time=%2.1f'%start)

	#plt.draw()
	im_index = 1
	plt.savefig("u_t_"+str(im_index)+".png")
	
	u_old = u0.copy()
	v_old = v0.copy()
	
	t=start
	
	# start loop
	for tout in np.arange(start+step,finish+step,step):
	    while t < tout:
			
			u_new = step_u(u_old, v_old, dt, dx2)
			v_new = step_v(u_old, v_old, dt, dx2)
			
			u_old = u_new.copy()
			v_old = v_new.copy()
						
			t+=dt
		
			
	    title.set_text('time=%.1f'%(t))
	    im.set_data((u_new.real).T)
	    #im.figure.canvas.draw()
	    im_index+=1
	    plt.savefig("u_t_"+str(im_index)+".png")
	    


def step_u(u_old, v_old, dt, dx2):
	"""
	Implement the step of equation \dot(u) = u - u^3 - v
	"""
	return u_old + dt*(u_old - u_old**3 - v_old + laplacian(u_old, dx2))

def step_v(u_old, v_old, dt, dx2):
	"""
	Implement the step of equation \dot(v) = e(u + a)
	"""
	return v_old + dt*(par['e']*(u_old + par['a1']*v_old + par['a0']) + par['delta']*laplacian(v_old,dx2))

def laplacian(u_old, h2):
	ui = u_old.copy()
	return (np.roll(ui,1,axis=0) + np.roll(ui,-1,axis=0) + np.roll(ui,1,axis=1) + np.roll(ui,-1,axis=1) -4*ui)/h2


if __name__ == "__main__":
	main()

