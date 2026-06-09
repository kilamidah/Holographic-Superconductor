# Holographic-Superconductor
This repository contains a high-precision numerical simulator for a holographic superconductor, developed in Python. Based on the AdS/CFT correspondence (gauge/gravity duality), the project models a (2+1)-dimensional strongly correlated superconducting system by solving the coupled nonlinear field equations of a (3+1)-dimensional classical gravitational dual. The model demonstrates the spontaneous breaking of a global U(1) symmetry via the formation of scalar "hair" around an Anti-de Sitter (AdS) planar black hole, signaling a second-order phase transition to a superconducting state when the temperature drops below a critical threshold, $T_c$
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
