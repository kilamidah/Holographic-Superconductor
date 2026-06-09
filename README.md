# Holographic Superconductor
This repository contains a high-precision numerical simulator for a holographic superconductor, developed in Python. Based on the AdS/CFT correspondence (gauge/gravity duality), the project models a (2+1)-dimensional strongly correlated superconducting system by solving the coupled nonlinear field equations of a (3+1)-dimensional classical gravitational dual. The model demonstrates the spontaneous breaking of a global U(1) symmetry via the formation of scalar "hair" around an Anti-de Sitter (AdS) planar black hole, signaling a second-order phase transition to a superconducting state when the temperature drops below a critical threshold, $T_c$.
# Theoretical Background
Conventional BCS theory struggles to describe strongly correlated, high-temperature superconductors (like cuprates). The AdS/CFT correspondence provides a powerful workaround by mapping strongly coupled quantum field theories at a boundary to weakly coupled classical gravity theories in a bulk spacetime with one extra dimension.

In this model:
 - The Temperature (T) of the superconductor is identified with the Hawking temperature of the black hole.
 - The Superconducting Condensate $(\langle O_2\rangle)$ corresponds to the normalizable mode of a charged scalar field in the bulk.
 - The Chemical Potential ($\mu$) and Charge Density ($\rho$) are extracted from the asymptotic behavior of the bulk gauge field.
# Key-Features
 - Numerical Pipeline: Converts complex boundary-value problems into efficient, solvable initial-value problems.
 - Stable Integration: Utilizes implicit Runge-Kutta (Radau) methods to handle stiffness near the black hole horizon.
 - Physical Validation: Reproduces established benchmark results, confirming accuracy in modeling spontaneous symmetry breaking.
# Technical Stack
Language = Python

Libraries = NumPy, SciPy, Matplotlib
# How To Run
Download the necessary python project file. This is named "Holographic Superconductor.py"

Before running make sure the three necessary libraries are all correctly installed. Run below codes individually.
 1. "pip install numpy"
 3. "pip install scipy"
 4. "pip install matplotlib"

Now run the main python file.
# Academic Context
This work was completed as part of my BSc in Theoretical Physics. It demonstrates the bridge between general relativity, condensed matter physics, and computational science.
