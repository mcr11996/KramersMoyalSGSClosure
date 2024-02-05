# -*- coding: utf-8 -*-
"""
Created on Mon May  8 20:39:39 2023

@author: Molly Ross
    This code is based off of the burgersLES code from Jeremy Gibbs' pyBurgers repository
    https://github.com/jeremygibbs/pyBurgers
    Additions have been made for a KM closure calculated from DNS
"""
import time
from sys import stdout
import numpy as np
import netCDF4 as nc
from burgers import Utils, Settings, BurgersLES
from KM_utils import find_KM_fit_coeffs, KM

utils = Utils()

def findTau(u_ss, delta_f=1,len_x=2*np.pi):
    """
    Find tau from DNS (u_ss) to put into Burgers Eq.

    Parameters
    ----------
    u_ss : TYPE numpy array
        Velocity series (spatial).
    delta_f : TYPE, integer
        Filter size ratio. The default is 1.
    len_x : TYPE, float
        Length of entire spatial domain. The default is 2*np.pi.

    Returns
    -------
    tau : TYPE numpy array
        Numpy array with length = len(u_ss)/delta_f.

    """
    uf = utils.filterDown(u_ss,delta_f)
    u2 = u_ss*u_ss
    term1 = utils.filterDown(u2,delta_f)
    term2 = uf*uf
    tau = term1 - term2
    dx = len_x/(len(tau))
    der = utils.derivative(tau,dx)
    return tau, der['dudx']


# LES solver
def main():

    # Let's time this thing
    t1 = time.time()

    # A nice welcome message
    print("##############################################################")
    print("#                                                            #")
    print("#                   Welcome to pyBurgers                     #")
    print("#      A fun tool to study turbulence using DNS and LES      #")
    print("#                                                            #")
    print("##############################################################")
    print("[pyBurgers: Info] \t You are running in LES mode")

    # Instantiate helper classes
    print("[pyBurgers: Setup] \t Reading input settings")
    utils    = Utils()
    settings = Settings('namelist.json')
    
    nxDNS = settings.nxDNS
    nxLES = 512
    mp    = int(nxLES/2)
    dx    = 2*np.pi/nxLES
    dt    = settings.dt
    nt    = settings.nt
    visc  = settings.visc
    damp  = settings.damp

    # Define the filename for the DNS training data
    dns_fname = 'pyBurgersDNS.nc'
    
    dns_data =  nc.Dataset(dns_fname,'r')
    
    u_dns = dns_data['u'][:]
    dns_data.close()

    # Define the time step for the DNS data training
    dt_DNS = 0.1
        
    # Calculate time series for tau from DNS to train KM model 
    for i in range(np.shape(u_dns)[0]):
        tau_i, dtaudx_i = findTau(u_dns[i,:],delta_f=int(nxDNS/nxLES))
        if 'tau_dns' in locals():
            tau_dns = np.vstack([tau_dns,tau_i])
            dtaudx_dns = np.vstack([dtaudx_dns,tau_i])
        if 'tau_dns' not in locals():
            tau_dns = tau_i
            dtaudx_dns = dtaudx_i
    
    # Find KM Coefficients
    bins, D1_e, D2_e = KM(tau_dns[:,80], lambda_1=False, dt=dt_DNS, num_bins=200)
    d1_coeffs, d2_coeffs = find_KM_fit_coeffs(bins,D1_e,D2_e)
    
    # Initiate random KM (This maintains the seed for the forcing function)
    eta = np.random.normal(0,1,[int(1000),int(nxLES)])
    
   
    # Initialize velocity field
    print("[pyBurgers: Setup] \t Initialzing velocity field")
    u = np.zeros(nxLES)

    # Initialize random number generator
    np.random.seed(1)

    # Place holder for right hand side
    rhsp = 0
    
    # Create output file
    # Some info commented out to speed up calculations and minimize output
    # file size.
    print("[pyBurgers: Setup] \t Creating output file")
    output = nc.Dataset('pyBurgersLES_KMfromDNS_Seed3.nc','w')
    output.description = "pyBurgers KM LES output"
    output.source = "M. Ross"
    output.history = "Created " + time.ctime(time.time())
    #output.setncattr("sgs","%d"%model)
    
    # Add dimensions
    output.createDimension('t')
    output.createDimension('x',nxLES)

    # Add variables
    out_t = output.createVariable("t", "f4", ("t"))
    out_t.long_name = "time"
    out_t.units = "s"
    out_x = output.createVariable("x", "f4", ("x"))
    out_x.long_name = "x-distance"
    out_x.units = "m"
    out_k = output.createVariable("tke", "f4", ("t"))
    out_k.long_name = "turbulence kinetic energy"
    out_k.units = "m2 s-2"
    # out_c = output.createVariable("C", "f4", ("t"))
    # out_c.long_name = "subgrid model coefficient"
    # out_c.units = "--"
    # out_ds = output.createVariable("diss_sgs", "f4", ("t"))
    # out_ds.long_name = "subgrid dissipation"
    # out_ds.units = "m2 s-3"
    # out_dm = output.createVariable("diss_mol", "f4", ("t"))
    # out_dm.long_name = "molecular dissipation"
    # out_dm.units = "m2 s-3"
    # out_ep = output.createVariable("ens_prod", "f4", ("t"))
    # out_ep.long_name = "enstrophy production"
    # out_ep.units = "s-3"
    # out_eds = output.createVariable("ens_diss_sgs", "f4", ("t"))
    # out_eds.long_name = "subgrid enstrophy dissipation"
    # out_eds.units = "s-3"
    # out_edm = output.createVariable("ens_diss_mol", "f4", ("t"))
    # out_edm.long_name = "molecular enstrophy dissipation"
    # out_edm.units = "s-3"
    out_u = output.createVariable("u", "f4", ("t","x"))
    out_u.long_name = "velocity"
    out_u.units = "m s-1"

    # Write x data
    out_x[:] = np.arange(0,2*np.pi,dx)
 
    # Time loop
    save_t = 0
    #dtaudx = np.zeros(nxLES)
    # Initiate tau value
    tau = tau_dns[0,:]
    tau_deriv = utils.derivative(tau,dx)
    dtaudx = tau_deriv['dudx']
    for t in range(int(nt)):
        
        # Update progress
        if (t==0 or (t+1)%1000==0):
            stdout.write("\r[pyBurgers: LES] \t Running for time %07d of %d"%(t+1,int(nt)))
            stdout.flush()
        
        # Compute derivatives
        derivs = utils.derivative(u,dx)
        dudx   = derivs['dudx']
        du2dx  = derivs['du2dx']
        d2udx2 = derivs['d2udx2']
        #d3udx3 = derivs['d3udx3']

        # Add fractional Brownian motion (FBM) noise
        fbm  = utils.noise(0.75,nxDNS)
        #fbmf = utils.noise(0.75,nxLES)
        fbmf = utils.filterDown(fbm,int(nxDNS/nxLES))

        # # compute subgrid terms from KM
        #tau = findTau(u, delta_f=1,len_x=2*np.pi)
        drift = (tau)*d1_coeffs[0]+d1_coeffs[1]
        diffusion = (tau)**2*d2_coeffs[0]+(tau)*d2_coeffs[1]+d2_coeffs[2]
        #diffusion = dtaudx**2*d2_coeffs[0]+dtaudx*d2_coeffs[1]+d2_coeffs[2]
        tau = tau + dt*drift + np.sqrt(2*diffusion*dt)*eta[int(np.remainder(t,1000)-0),:]#np.random.normal(0,1,nxLES)
        tau_deriv = utils.derivative(tau,dx)
        dtaudx = tau_deriv['dudx']

        # Compute right hand side
        rhs = visc * d2udx2 - 0.5*du2dx + np.sqrt(2*damp/dt)*fbmf - 0.5*dtaudx
        
        # Time integration
        if t == 0:
            # Euler for first time step
            u_new = u + dt*rhs
        else:
            # 2nd-order Adams-Bashforth
            u_new = u + dt*(1.5*rhs - 0.5*rhsp)
        
        # Set Nyquist to zero
        fu_new     = np.fft.fft(u_new)
        fu_new[mp] = 0
        u_new      = np.real(np.fft.ifft(fu_new))
        u          = u_new
        rhsp       = rhs

        # Output to file every 1000 time steps (0.1 seconds)
        if ((t+1)%1000==0):
            
            # Fix the tau value to the DNS data every available time step
            tau = tau_dns[save_t,:]
            
            # Kinetic energy
            tke  = 0.5*np.var(u)

            # Dissipation
            #diss_sgs = np.mean(-tau*dudx)
            #diss_mol = np.mean(visc*dudx**2)

            # Enstrophy
            #ens_prod = np.mean(dudx**3)
            #ens_dsgs = np.mean(-tau*d3udx3)
            #ens_dmol = np.mean(visc*d2udx2**2)
            
            # Save to disk
            out_t[save_t]   = (t+1)*dt 
            out_k[save_t]   = tke
            #out_c[save_t]   = coeff
            #out_ds[save_t]  = diss_sgs
            #out_dm[save_t]  = diss_mol
            #out_ep[save_t]  = ens_prod
            #out_eds[save_t] = ens_dsgs
            #out_edm[save_t] = ens_dmol
            out_u[save_t,:] = u
            save_t += 1
            #u = u_dns[save_t,::int(nxDNS/nxLES)]
    
    # Time info
    t2 = time.time()
    tt = t2 - t1
    print("\n[pyBurgers: LES] \t Done! Completed in %0.2f seconds"%tt)
    print("##############################################################")

if __name__ == "__main__":
    main()
