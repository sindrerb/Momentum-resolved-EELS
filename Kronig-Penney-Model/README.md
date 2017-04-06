# Kronig-Penney-Model
This program generates a bandstructure with corresponding wavefunctions from a lattice with delta potentials. 
The program is heavily inspired by the work of Robert B. Laughlin (http://large.stanford.edu/courses/2007/ap272/laughlin2/)

The real and reciprocal lattice are defined by a CELLFILE. A plane wave basis is constructed by the reciprocal lattice vectors in combination with a energy cut off  and saved to a BASISFILE.

Initially the program assumes an empty lattice thus nearly free electrons. The system is perturbed with deltapotentials at each lattice point, the eigenenergies and eigenstates of the perturbed system is written to a WAVEFILE.

The WAVEFILE may be used as an initial state for a new perturbation, allowing to disturb the system until a desired band structure is achieved.

