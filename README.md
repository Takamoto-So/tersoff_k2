pair tersoff_k2
==============

_So Takamoto_

LAMMPS implementation of the interatomic potential for Si-C system.  

The related works have published in Physical Review B and Journal of Applied Physics:
So Takamoto, Takahiro Yamasaki, Jun Nara, Takahisa Ohno, Chioko Kaneta, Asuka Hatano, and Satoshi Izumi, "Atomistic mechanism of graphene growth on SiC substrate: Large-scale molecular dynamics simulation based on a new charge-transfer bond-order type potential", Physical Review B, 97, 125411 (2018).
DOI: <https://doi.org/10.1103/PhysRevB.97.125411>
ArXiv: <https://arxiv.org/abs/1802.07871>
So Takamoto, Takahiro Yamasaki, Takahisa Ohno, Chioko Kaneta, Asuka Hatano, and Satoshi Izumi, "Elucidation of the atomic-scale mechanism of the anisotropic oxidation rate of 4H-SiC between the (0001) Si-face and (000-1) C-face by using a new Si-O-C interatomic potential", Journal of Applied Physics.

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

