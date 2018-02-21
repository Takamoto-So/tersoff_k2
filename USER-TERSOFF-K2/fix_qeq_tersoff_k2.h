/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author:
   Hasan Metin Aktulga, Purdue University - fix qeq/reax implementation
    (now at Lawrence Berkeley National Laboratory, hmaktulga@lbl.gov)
   So Takamoto, The University of Tokyo - fix qeq/tersoff_k(2) implementation, which is based on fix qeq/reax

   Please cite the related publication:
     fix qeq/reax:
     H. M. Aktulga, J. C. Fogarty, S. A. Pandit, A. Y. Grama,
     "Parallel Reactive Molecular Dynamics: Numerical Methods and
     Algorithmic Techniques", Parallel Computing, in press.

     fix qeq/tersoff_k(2):
     S. Takamoto, T. Kumagai, T. Yamasaki, T. Ohno, C. Kaneta, A. Hatano, S. Izumi,
     "Charge-transfer interatomic potential for investigation of the thermal-oxidation growth process of silicon", J. of Appl. Phys.

------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(qeq/tersoff/k2,FixQEqTersoff_k2)

#else

#ifndef LMP_FIX_QEQ_TERSOFF_K2_H
#define LMP_FIX_QEQ_TERSOFF_K2_H

#include "fix.h"
#include "pair_tersoff_k2.h"

namespace LAMMPS_NS {

class FixQEqTersoff_k2 : public Fix {
 public:
  FixQEqTersoff_k2(class LAMMPS *, int, char **);
  ~FixQEqTersoff_k2();
  int setmask();
  void init();
  void init_list(int,class NeighList *);
  void init_storage();
  void setup_pre_force(int);
  void pre_force(int);

  void setup_pre_force_respa(int, int);
  void pre_force_respa(int, int, int);

  void min_setup_pre_force(int);
  void min_pre_force(int);

  int matvecs;
  double qeq_time;

 private:
  int nevery,/*reaxflag,*/ tersoff_kflag;
  //int n, N, m_fill;
  int n_cap, nmax, m_cap;
  int pack_flag;
  int nlevels_respa;
  class NeighList *list;
  //class PairReaxC *reaxc;
  class PairTersoff_k2 *tersoff_k;

  double swa, swb;      // lower/upper Taper cutoff radius
  double Tap[8];        // Taper function
  double tolerance;     // tolerance for the norm of the rel residual in CG

  double *chi,*eta,*gamma;  // qeq parameters
  std::vector<double> chi_mod, eta_mod;
  //double *chi_mod, *eta_mod;
  double **shld;

  // fictitious charges

  //double *s, *t;
  //double **s_hist, **t_hist;
  std::vector<double> s, t;
  std::vector<double> s_hist, t_hist;
  const int nprev;

  TERSOFF_K2::sparse_matrix H;
  //double *Hdia_inv;
  //double *b_s, *b_t;
  //double *b_prc, *b_prm;
  std::vector<double> Hdia_inv;
  std::vector<double> b_s, b_t;
  std::vector<double> b_prc, b_prm;

  //CG storage
  //double *p, *q, *r, *d;
  std::vector<double> p_vec, q_vec, r_vec, d_vec;

  //GMRES storage
  //double *g,*y;
  //double **v;
  //double **h;
  //double *hc, *hs;

  void pertype_parameters(char*);
  void init_shielding();
  void init_taper();
  void allocate_storage();
  void deallocate_storage();
  void reallocate_storage();

  void copy_pairdata();

  void init_matvec();
  void calculate_Q();

  int CG(const TERSOFF_K2::sparse_matrix&,double*,double*);
  //int GMRES(double*,double*);
  void sparse_matvec(const TERSOFF_K2::sparse_matrix&,double*,double*);

  int pack_forward_comm(int, int *, double *, int, int *);
  void unpack_forward_comm(int, int, double *);
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
  double memory_usage();
  void grow_arrays(int);
  void copy_arrays(int, int, int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);

  double parallel_norm( double*, int );
  double parallel_dot( double*, double*, int );
  double parallel_vector_acc( double*, int );

  double norm(double*,int);
  void vector_sum(double*,double,double*,double,double*,int);
  void vector_scale(double*,double,double*,int);
  double dot(double*,double*,int);
  void vector_add(double*, double, double*,int);
};

}

#endif
#endif
