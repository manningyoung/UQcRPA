# UQcRPA

A MATLAB package for computing the screened Coulomb interaction (Hubbard _U_) in solids via the constrained random phase approximation (cRPA). Interfaces with [BerkeleyGW](http://www.berkeleygw.org/) and [wannier90](http://www.wannier.org/).

**(!) Very early release (!)**  
Contact the author for usage instructions.

---

## Input data description

UQcRPA requires the following input data from **BerkeleyGW**:
* `chi0mat.h5`, `chimat.h5`: Full polarisability matrix.
* `chi0mat_a.h5`, `chimat_a.h5`: Polarisability matrix for the active subspace.
* `eps0mat.h5`, `epsmat.h5`: Full inverse dielectric matrix (required for bare Coulomb interaction).
* `epsilon.log`: Logfile from epsilon calculation (required for reciprocal lattice vectors).

As well as the following from **wannier90**:
* `seedname_u.mat`: Unitary MLWF transformation matrix.
* Set of `UNKp.s`: Periodic part of the Bloch states used for Wannier construction.

## Tips

1. Save the transformed Bloch states `unk` as `unk.mat` in the data directory to speed up repeat calculations.