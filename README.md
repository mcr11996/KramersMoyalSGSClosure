# Kramers Moyal Burgulence Closure
This repository contains code to provide a data-driven closure for the subgrid scale term in Burgers turbulence (Burgulence) using the Kramers-Moyal expansion method.

## Stochastic Burgers Equation

This code is built off of [Jeremy Gibbs' pyBurgers repository](https://github.com/jeremygibbs/pyBurgers). A more detailed description of the numerical method for solving the stochastic Burgers equation for both DNS and LES can be found there.

## Kramers-Moyal Expansion Method

The Kramers-Moyal (KM) expansion method is a statistical mechanics-based method for describing the evolution of a pdf in time.

The generic equation for some state variable, $x$, is given by

$\frac{dt}{dx}=D_1(x,t)+\sqrt{D_2(x,t)}dW$

This repository assumes that the 

# Using the Code

## Code dependencies

This code is written in Python 3 and is compatible with Windows, Linux, and Mac operating systems.

Required libraries include . To install,

```
pip install scikit.sklearn
```

# Licensing

This code is licensed under GNU GPL-3.0.

# Contribute

Suggestions to improve this repository are always welcome.

Errors can be reported to

# Citation

If using this repository, please cite


