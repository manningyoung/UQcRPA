# UQcRPA

A MATLAB package for computing the screened Coulomb interaction (Hubbard _U_) in solids via the constrained random phase approximation (cRPA). Interfaces with [BerkeleyGW](http://www.berkeleygw.org/) and [wannier90](http://www.wannier.org/).

**Version 0.1** (May 2017)  
M. C. Young, A. C. Jacko, B. J. Powell.  
Tested on MATLAB R2017a with BerkeleyGW 1.2.0 and wannier90 2.1.0.  
Known bugs - use with caution.

## Quick start guide
The package is run from the `main` script. In this file, specify the following:
* `datadir`: Location of the input data (described below). Defaults to the UQcRPA directory.
* `wann_bands`: Band indices corresponding to the active subspace (Wannier bands).
* `iband`, `jband`: The matrix element U_{iband, jband} to be calculated.


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

* Save the transformed Bloch states `unk` as `unk.mat` in the data directory to speed up repeat calculations.