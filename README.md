pair tersoff_k2
==============

_So Takamoto_

LAMMPS implementation of the interatomic potential for Si-C system.  

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

Please understand that we cannot answer questions.
