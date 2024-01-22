# Kramers Moyal Burgulence Closure
This repository contains code to provide a data-driven closure for the subgrid scale term in Burgers turbulence (Burgulence) using the Kramers-Moyal expansion method. All relevant code can be found in the Code folder.

## Stochastic Burgers Equation

This code is built off of [Jeremy Gibbs' pyBurgers repository](https://github.com/jeremygibbs/pyBurgers). A more detailed description of the numerical method for solving the stochastic Burgers equation for both DNS and LES can be found there.

The equation for solving Burgers equation using LES is defined by

$\frac{\partial \bar{u}}{\partial t}+\bar{u}\frac{\partial \bar{u}}{\partial x} = \nu \frac{\partial^2\bar{u}}{\partial x^2} + \bar{f}(x,t)-\frac{1}{2}\frac{\partial \tau}{\partial x}$

The objective of this code is to generate a model for $\tau$ in the subgrid-scale term.

## Kramers-Moyal Expansion Method

The Kramers-Moyal (KM) expansion method is a statistical mechanics-based method for describing the evolution of a pdf in time.

The generic equation for some state variable, $x$, is given by

$\frac{dt}{dx}=D_1(x,t)+\sqrt{D_2(x,t)}dW$

The value for the subgrid-scale (SGS) term, $\tau$, can be calculated from the resultant DNS solution.

$\tau = \overline{uu} - \bar{u}\bar{u}$

The values for $D_1$ and $D_2$ can be calculated from time series of $\tau$ calculated from DNS.

The value for $\tau$ can be found on the fly by solving the SDE numerically using the Euler-Marayuma method.

$\tau_{n+1} = \tau_n + D_1(\tau_n)dt
    + D_2(\tau_n)dW_n+\frac{1}{2}(D_2D_2')(\tau_n)(dW_n^2-dt)$

See the [the cited paper](https://www.researchgate.net/profile/Hitesh-Bindra/publication/376766325_Development_of_a_subgrid-scale_model_for_Burgers_turbulence_using_statistical_mechanics-based_methods/links/65999ab53c472d2e8eb968a9/Development-of-a-subgrid-scale-model-for-Burgers-turbulence-using-statistical-mechanics-based-methods.pdf) for a more thorough derivation of the math involved.

# Using the Code

This is a data-driven closure for Burgers turbulence. To obtain the initial data, run pyBurgersDNS.py. This will produce a netCDF4 (.nc) file with the results.

Next, run BurgersLESKMfromDNS.py. This will produce an LES 

## Code dependencies

This code is written in Python 3 and is compatible with Windows, Linux, and Mac operating systems.

Required libraries include numpy, scipy, scikit.sklearn, and netCDF4. To install,

```
pip install scikit.sklearn
```

and

```
pip install netCDF4
```

# Licensing

This code is licensed under GNU GPL-3.0.

# Contribute

Suggestions to improve this repository are always welcome. Please report any errors through Github.

# Citation

If using this repository, please cite

Ross, Molly, and Hitesh Bindra. "Development of a subgrid-scale model for Burgers turbulence using statistical mechanics-based methods." Physics of Fluids 35.12 (2023).


