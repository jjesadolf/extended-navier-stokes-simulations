# extended-navier-stokes-simulations
MATLAB numerical simulations exploring nonlinear thermo-vortical dynamics, stability regimes and extended Navier–Stokes inspired models.

# Trying to Understand How Nonlinear Feedback Shapes Fluid Dynamics

## Numerical Exploration of Thermo-Vortical Dynamics Inspired by Extended Navier–Stokes Frameworks

### Author

Jesus Adolfo Rossi

---

# Abstract

This project presents an exploratory numerical framework developed in MATLAB to investigate the interaction between nonlinear feedback, vorticity evolution, energy redistribution and stability regimes in simplified thermo-vortical systems inspired by extended Navier–Stokes dynamics.

The objective of the work is not to propose a definitive fluid model, but rather to study how additional nonlinear thermo-vortical coupling mechanisms influence the emergence of bounded oscillatory regimes, spectral energy distributions and reduced phase-space structures.

Several numerical configurations were investigated, including:

* classical Navier–Stokes inspired dynamics,
* thermo-vortical feedback,
* memory-driven interactions,
* anisotropic coupling structures,
* nonlinear saturation mechanisms,
* local stability diagnostics.

The simulations reveal that relatively small modifications in nonlinear feedback terms can generate significantly different dynamical behaviors, ranging from bounded spiral-like trajectories to stronger oscillatory and destabilizing regimes.

The project combines:

* numerical simulation,
* reduced-order nonlinear dynamics,
* spectral analysis,
* effective vorticity diagnostics,
* phase-space reconstruction,
* stability-oriented interpretation.

This repository is intended as an open exploratory framework for studying nonlinear dynamical behaviors emerging from thermo-vortical coupling.

---

# 1. Motivation

One of the most fascinating aspects of fluid dynamics is the emergence of highly complex structures from relatively compact mathematical descriptions.

As nonlinear interactions become dominant, the evolution of the system may develop:

* coherent vortical structures,
* oscillatory regimes,
* intermittent amplification,
* multiscale energy redistribution,
* localized instability.

This project was born from the desire to better understand how additional nonlinear thermo-vortical feedback mechanisms may influence these behaviors.

Rather than searching for a final or exact physical theory, the objective is exploratory:

to observe how the interaction between energy, vorticity, memory and nonlinear coupling reorganizes the global and local dynamics of the system.

Particular attention was given to:

* bounded vs destabilizing regimes,
* oscillatory structures,
* phase-space evolution,
* spectral signatures,
* local stability indicators.

---

# 2. Conceptual Framework

The simulations are inspired by incompressible Navier–Stokes-type dynamics extended through additional thermo-vortical coupling structures.

The conceptual idea behind the framework is that local vorticity evolution may interact with:

* thermal gradients,
* nonlinear feedback,
* memory effects,
* anisotropic structures,
* saturation mechanisms.

A representative thermo-vortical contribution used throughout the exploratory analyses is:

[
F_{\eta}(T,\omega)=-\eta\frac{\nabla T\times\omega}{1+c_{\omega}|\omega|^2}
]

where:

* (\omega) represents vorticity,
* (T) represents a thermal-like field,
* (\eta) controls coupling intensity,
* (c_{\omega}) introduces nonlinear saturation.

The saturation denominator was introduced to prevent unrealistic divergence and to investigate whether bounded nonlinear regimes could emerge spontaneously.

---

# 3. Numerical Objectives

The numerical analyses were designed to investigate:

## 3.1 Energy Redistribution

How nonlinear feedback structures redistribute energy across time and dynamical modes.

Questions explored:

* Do some configurations dissipate energy more efficiently?
* Which regimes generate persistent oscillatory exchange?
* Can nonlinear coupling confine energetic growth?

---

## 3.2 Effective Vorticity Evolution

Vorticity was treated as a central dynamical quantity because coherent structures in many fluid systems emerge through vortical interactions.

The analyses investigated:

* growth,
* saturation,
* oscillatory behavior,
* boundedness,
* local amplification.

---

## 3.3 Reduced Phase-Space Dynamics

Reduced phase-space trajectories were reconstructed to visualize the geometry of the nonlinear evolution.

The objective was to observe whether different thermo-vortical structures generate:

* compact bounded trajectories,
* spiral-like regimes,
* expanding oscillatory dynamics,
* attractor-like structures,
* destabilizing transitions.

---

## 3.4 Spectral Analysis

Fourier-based analyses were introduced to study how different regimes redistribute energy across frequencies.

Particular attention was given to:

* low-frequency concentration,
* spectral broadening,
* multiscale activity,
* oscillatory signatures.

---

## 3.5 Local Stability Diagnostics

Additional diagnostics were progressively introduced to distinguish physically meaningful bounded structures from purely numerical amplification.

These included:

* effective stability indicators,
* thermo-vortical control ratios,
* localized high-vorticity analysis,
* comparison between classical vortex stretching and thermo-vortical contributions.

---

# 4. Phase-Space Trajectories

One of the most visually interesting aspects of the simulations was the emergence of qualitatively different phase-space geometries.

Some trajectories remained compact and spiral-like, suggesting bounded oscillatory regimes.

Others expanded into wider structures associated with stronger nonlinear coupling and more persistent oscillatory activity.

These behaviors were particularly interesting because they were not explicitly programmed.

They emerged spontaneously from the nonlinear interaction between:

* feedback,
* vorticity evolution,
* memory,
* anisotropy,
* saturation mechanisms.

The reduced trajectories suggest that even relatively small modifications in nonlinear coupling may reorganize the global dynamics in substantial ways.

---

# 5. Energy Evolution

The energy evolution plots revealed important qualitative differences between the tested configurations.

Some regimes converged toward nearly stationary energetic states.

Others maintained persistent oscillatory exchanges between nonlinear modes.

In stronger coupling regimes, the redistribution of energy became increasingly dynamic and less confined.

Anisotropic and memory-driven configurations frequently generated richer oscillatory behavior compared to simpler baseline models.

These results suggest that nonlinear thermo-vortical feedback may significantly influence:

* confinement,
* dissipation,
* energy transfer mechanisms,
* long-term oscillatory evolution.

---

# 6. Effective Vorticity Dynamics

The effective vorticity diagnostics became one of the central tools of the project.

Monitoring the evolution of vorticity allowed the identification of:

* stabilization,
* saturation,
* oscillatory amplification,
* localized nonlinear growth.

Some configurations remained strongly bounded over long timescales.

Others exhibited increasingly wide oscillatory behavior associated with stronger nonlinear coupling.

These analyses highlighted how subtle changes in feedback structure may substantially alter the balance between confinement and amplification.

---

# 7. Spectral Analysis

The spectral analyses revealed that different configurations develop distinct frequency-domain signatures.

Some regimes concentrated most of the energy into low-frequency modes.

Others generated broader spectra associated with more complex multiscale activity.

The spectral broadening observed in stronger nonlinear regimes suggests that the thermo-vortical interactions may redistribute activity across multiple scales rather than maintaining simple coherent oscillatory behavior.

This aspect may be particularly relevant for future investigations involving:

* turbulence,
* nonlinear transport,
* multiscale coupling,
* reduced-order fluid dynamics.

---

# 8. Local Stability Diagnostics

To better understand the local organization of the dynamics, additional stability-oriented diagnostics were introduced.

These included:

* effective stability indicators,
* critical-region analysis,
* thermo-vortical control ratios,
* comparison between classical vortex stretching and thermo-vortical contributions.

One particularly interesting observation was that, in some regions, the thermo-vortical contribution appeared to partially counteract classical stretching mechanisms.

This does not constitute a rigorous physical conclusion, but it suggests that nonlinear thermo-vortical feedback may reorganize local amplification processes in nontrivial ways.

The purpose of these diagnostics was to move beyond purely global observables and investigate where destabilization or confinement was actually emerging inside the evolving nonlinear dynamics.

---

# 9. Interpretation and Limitations

The project remains exploratory and qualitative in nature.

The simulations should not be interpreted as rigorous proofs of new fluid theories.

Instead, they should be viewed as a numerical framework designed to investigate how nonlinear feedback interactions generate different dynamical structures.

Several important limitations remain:

* reduced-order simplifications,
* numerical sensitivity,
* dependence on parameter selection,
* absence of rigorous convergence analysis,
* lack of formal physical validation.

Future improvements may include:

* Lyapunov exponent analysis,
* Poincaré sections,
* convergence studies,
* comparison with known nonlinear attractors,
* improved stability diagnostics,
* higher-resolution numerical schemes.

---

# 10. Future Directions

Potential future developments include:

* pseudo-spectral 3D simulations,
* stronger PDE coupling,
* reduced-order attractor analysis,
* anisotropic geometric structures,
* turbulence-oriented diagnostics,
* nonlinear stability characterization,
* data-driven analysis of oscillatory regimes.

Another important direction would be the comparison between the observed structures and known nonlinear dynamical behaviors in:

* turbulence transitions,
* bounded attractors,
* multiscale oscillatory systems,
* nonlinear transport phenomena.

---

# 11. Final Thoughts

Perhaps the most fascinating aspect of these simulations is that many structures were never explicitly designed.

They emerge spontaneously from the nonlinear interaction between energy redistribution, vorticity evolution, memory and feedback.

At times, numerical modeling feels less like forcing equations to behave and more like observing a dynamical system reveal its own internal tendencies.

That curiosity is ultimately what motivated this project.

---

# Tools and Technologies

* MATLAB
* Numerical Simulation
* Fluid Dynamics
* Nonlinear Dynamics
* Stability Analysis
* Spectral Analysis
* Reduced-Order Modeling

---

# Repository Notes

This repository is continuously evolving.

Feedback, critiques and suggestions are always welcome, particularly from researchers and engineers working in:

* CFD,
* turbulence,
* nonlinear systems,
* applied mathematics,
* numerical simulation,
* computational physics.
