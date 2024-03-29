#/usr/bin/env python
# THIS FILE OBTAINED FROM JEREMY GIBBS' PYBURGERS PROJECT
# https://github.com/jeremygibbs/pyBurgers
import time
from sys import stdout
import numpy as np
import netCDF4 as nc
from burgers import Utils, Settings

# DNS run loop
def main():

    # let's time this thing
    t1 = time.time()

    # a nice welcome message
    print("##############################################################")
    print("#                                                            #")
    print("#                   Welcome to pyBurgers                     #")
    print("#      A fun tool to study turbulence using DNS and LES      #")
    print("#                                                            #")
    print("##############################################################")
    print("[pyBurgers: Info] \t You are running in DNS mode")

    # instantiate helper classes
    print("[pyBurgers: Setup] \t Reading input settings")
    utils    = Utils()
    settings = Settings('namelist.json')

    # input settings
    nx   = settings.nxDNS
    mp   = int(nx/2)
    dx   = 2*np.pi/nx
    dt   = settings.dt
    nt   = settings.nt
    visc = settings.visc
    damp = settings.damp
    
    # initialize velocity field
    print("[pyBurgers: Setup] \t Initialzing velocity field")
    u = np.zeros(nx)

    # initialize random number generator
    np.random.seed(1)

    # place holder for right hand side
    rhsp = 0
  
    # create output file
    print("[pyBurgers: Setup] \t Creating output file")
    output = nc.Dataset('pyBurgersDNS.nc','w')
    output.description = "pyBurgers DNS output"
    output.source = "Jeremy A. Gibbs"
    output.history = "Created " + time.ctime(time.time())
    
    # add dimensions
    output.createDimension('t')
    output.createDimension('x',nx)

    # add variables
    out_t = output.createVariable("t", "f4", ("t"))
    out_t.long_name = "time"
    out_t.units = "s"
    out_x = output.createVariable("x", "f4", ("x"))
    out_x.long_name = "x-distance"
    out_x.units = "m"
    out_k = output.createVariable("tke", "f4", ("t"))
    out_k.long_name = "turbulence kinetic energy"
    out_k.units = "m2 s-2"
    out_u = output.createVariable("u", "f4", ("t","x"))
    out_u.long_name = "velocity"
    out_u.units = "m s-1"

    # write x data
    out_x[:] = np.arange(0,2*np.pi,dx)

    # time loop
    save_t = 0
    for t in range(int(nt)):

        # update progress
        if (t==0 or (t+1)%1000==0):
            stdout.write("\r[pyBurgers: DNS] \t Running for time %07d of %d"%(t+1,int(nt)))
            stdout.flush()
        
        # compute derivatives
        derivs = utils.derivative(u,dx)
        du2dx  = derivs['du2dx']
        d2udx2 = derivs['d2udx2'] 

        # add fractional Brownian motion (FBM) noise
        fbm = utils.noise(0.75,nx)

        # compute right hand side
        rhs = visc * d2udx2 - 0.5*du2dx + np.sqrt(2*damp/dt)*fbm
        
        # time integration
        if t == 0:
            # Euler for first time step
            u_new = u + dt*rhs
        else:
            # 2nd-order Adams-Bashforth
            u_new = u + dt*(1.5*rhs - 0.5*rhsp)
        
        # set Nyquist to zero
        fu_new     = np.fft.fft(u_new)
        fu_new[mp] = 0
        u_new      = np.real(np.fft.ifft(fu_new))
        u          = u_new
        rhsp       = rhs

        # output to file every 1000 time steps (0.1 seconds)
        if ((t+1)%1000==0):          
            
            # kinetic energy
            tke  = 0.5*np.var(u)
            
            # save to disk
            out_t[save_t]   = (t+1)*dt 
            out_k[save_t]   = tke
            out_u[save_t,:] = u
            save_t += 1

    # time info
    t2 = time.time()
    tt = t2 - t1
    print("\n[pyBurgers: DNS] \t Done! Completed in %0.2f seconds"%tt)
    print("##############################################################")


if __name__ == "__main__":
    main()
