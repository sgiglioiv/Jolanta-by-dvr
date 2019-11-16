# Jolanta-by-dvr
L^2 data analysis for resonances

Playground for different methods to find the complex energy of a resonance state (Eres = Er - i Gamma/2) from a series of calculations using square-integrable (L^2) wavefunctions.

I can think of four methods:
- complex scaling
- complex absorbing potebntials
- stabilization method (Hazi-Taylor--not extrapolation to zero stabilization, that is poor man's ACCC)
- ACCC or RAC = complex extrapolation of the lowest root

All four separate again along finer distinctions, but let's not go there.

For all four, the Hamiltonian is parametrized, and in step 1 of the calculation, the low-energy spectrum is repeatedly computed for many values of the parameter. The ACCC/RAC method represents an exception, as it only requires the lowest state, but the tradeoff is that only "low-energy" resonances can be found. Step 1 then yields a data set: energy_j(parameter_value_k).

In step 2, the data set is analysed to yield Eres, and the four methods differ in both the way the Hamiltionian is parametrized in step 1 and in the way the resulting data set is analyzed in step 2.

This project aims at comparing the four different approches:

- The Jolanta model potential is used
  1D:  V(x) = (a*x^2 - b) * exp(-c*x^2)
  3D:  V(r) = (a*r^2 - b) * exp(-c*r^2) + 0.5*l*(l+1)/r**2

- A simple sine DVR of the Hamiltonian is used

- Each step goes into a Jupyter notebook
