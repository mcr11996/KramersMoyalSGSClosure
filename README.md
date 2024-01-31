[![Contributors][contributors-shield]][contributors-url]
[![Forks][forks-shield]][forks-url]
[![Stargazers][stars-shield]][stars-url]
[![Issues][issues-shield]][issues-url]
[![GNU GPL-3.0 License][license-shield]][license-url]
[![Citation][citation-shield]][citation-url]


# Kramers Moyal Subgrid-scale Closure

This repository contains code to provide a data-driven closure for the subgrid scale term in Large Eddy Simulation (LES) of the momentum equation using the Kramers-Moyal expansion method. The current code is implemented in the 1D Burgers equation (Burger's turbulence). All relevant code can be found in the Code folder.

<!-- TABLE OF CONTENTS -->
# Table of Contents

<ol>
    <li>
      <a href="#stochastic-burgers-equation">Stochastic Burgers Equation</a>
      <ul>
        <li><a href="#built-with">Built With</a></li>
      </ul>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#roadmap">Roadmap</a></li>
    <li><a href="#contribute">Contribute</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#citation">Citation</a></li>
</ol>

<!--STOCHASTIC BURGERS EQUATION-->
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

<!--DEPENDENCIES-->
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


<!-- LICENSE -->
## License

Distributed under the GNU GPL-3.0 License. See `LICENSE.txt` for more information.

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- CONTRIBUTE -->
## Contribute

Contributions to this repository are always welcome. Any contributions you make are **greatly appreciated**.

If you have a suggestion that would make this project better, please fork the repo and create a pull request. You can also simply open an issue with the tag "enhancement".

1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature/AmazingFeature`)
3. Commit your Changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the Branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!--CITATION-->
# Citation

If using this repository, please cite

Ross, Molly, and Hitesh Bindra. "Development of a subgrid-scale model for Burgers turbulence using statistical mechanics-based methods." Physics of Fluids 35.12 (2023).

```
@article{ross2023development,
  title={Development of a subgrid-scale model for Burgers turbulence using statistical mechanics-based methods},
  author={Ross, Molly and Bindra, Hitesh},
  journal={Physics of Fluids},
  volume={35},
  number={12},
  year={2023},
  publisher={AIP Publishing}
}
```

[contributors-shield]: https://img.shields.io/github/contributors/mcr11996/KramersMoyalBurgulenceClosure?style=for-the-badge
[contributors-url]: https://github.com/mcr11996/KramersMoyalSGSClosure/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/mcr11996/KramersMoyalBurgulenceClosure?style=for-the-badge
[forks-url]: https://github.com/mcr11996/KramersMoyalBurgulenceClosure/network/members
[stars-shield]: https://img.shields.io/github/last-commit/mcr11996/KramersMoyalSGSClosure?style=for-the-badge
[stars-url]: https://github.com/mcr11996/KramersMoyalBurgulenceClosure
[issues-shield]: https://img.shields.io/github/issues/mcr11996/KramersMoyalSGSClosure?style=for-the-badge&color=red
[issues-url]: https://github.com/mcr11996/KramersMoyalSGSClosure/issues
[license-shield]: https://img.shields.io/badge/GPL-3?style=for-the-badge&label=LICENSE&color=green
[license-url]: https://github.com/mcr11996/KramersMoyalSGS/blob/master/LICENSE.txt
[citation-shield]: https://img.shields.io/badge/PUBLICATION-3?style=for-the-badge&logo=googlescholar&labelColor=grey&color=grey
[citation-url]: https://pubs.aip.org/aip/pof/article/35/12/125144/2930728
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555
[linkedin-url]: https://linkedin.com/in/molly-ross-48300186/
[product-screenshot]: images/screenshot.png
[Next.js]: https://img.shields.io/badge/next.js-000000?style=for-the-badge&logo=nextdotjs&logoColor=white
[Next-url]: https://nextjs.org/
[React.js]: https://img.shields.io/badge/React-20232A?style=for-the-badge&logo=react&logoColor=61DAFB
[React-url]: https://reactjs.org/
[Vue.js]: https://img.shields.io/badge/Vue.js-35495E?style=for-the-badge&logo=vuedotjs&logoColor=4FC08D
[Vue-url]: https://vuejs.org/
[Angular.io]: https://img.shields.io/badge/Angular-DD0031?style=for-the-badge&logo=angular&logoColor=white
[Angular-url]: https://angular.io/
[Svelte.dev]: https://img.shields.io/badge/Svelte-4A4A55?style=for-the-badge&logo=svelte&logoColor=FF3E00
[Svelte-url]: https://svelte.dev/
[Laravel.com]: https://img.shields.io/badge/Laravel-FF2D20?style=for-the-badge&logo=laravel&logoColor=white
[Laravel-url]: https://laravel.com
[Bootstrap.com]: https://img.shields.io/badge/Bootstrap-563D7C?style=for-the-badge&logo=bootstrap&logoColor=white
[Bootstrap-url]: https://getbootstrap.com
[JQuery.com]: https://img.shields.io/badge/jQuery-0769AD?style=for-the-badge&logo=jquery&logoColor=white
[JQuery-url]: https://jquery.com 


