pair tersoff_k2
==============

_So Takamoto_

LAMMPS implementation of the interatomic potential for Si-C system.  

The related work has accepted in Physical Review B:
Takamoto S. et al. (2018) Atomistic mechanism of graphene growth on SiC substrate: Large-scale molecular dynamics simulation based on a new charge-transfer bond-order type potential Physical Review B, American Physical Society. Forthcoming 2018.
ArXiv: <https://arxiv.org/abs/1802.07871>

An example of molecular dynamics simulation using this interatomic potential:
[![Graphene Growth Simulation (thermal decomposition of SiC)](http://img.youtube.com/vi/s5T1AEZ5G_0/0.jpg)](http://www.youtube.com/watch?v=s5T1AEZ5G_0)


Installation
------------

The `pair_style tersoff/k2`, the `fix qeq/tersoff/k2` are included
in the LAMMPS distribution as the USER-TERSOFF-K2 package.

To compile:

    cp -r 'USER-TERSOFF-K2' 'lammps/src'
    cd 'lammps/src'
    make yes-USER-TERSOFF-K2
    make 'machine'


Documentation
-------------

The usage of `pair_style` and `pair_coeff` is same to the original tersoff and tersoff/k potential.
For tersoff/k, see <https://github.com/Takamoto-So/tersoff_k>

The usage of `fix qeq/tersoff/k2` is same to the `fix qeq/tersoff/k` command.  

Other
-----

We have tested this package in LAMMPS 16-Feb-16 version.  
It may not work in some environments/versions.

