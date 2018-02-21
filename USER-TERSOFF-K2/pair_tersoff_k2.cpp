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
   Contributing author: Aidan Thompson (SNL) - original Tersoff implementation
                        So Takamoto (The University of Tokyo) - Hybrid Tersoff implementation (this file)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_tersoff_k2.h"
#include "atom.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "force.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

#include "math_const.h"

using namespace LAMMPS_NS;
using namespace LAMMPS_NS::TERSOFF_K2;
using namespace MathConst;

const int MAXLINE = 1024;
const int DELTA = 4;

/* ---------------------------------------------------------------------- */

PairTersoff_k2::PairTersoff_k2(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
  restartinfo = 0;
  one_coeff = 1;
  manybody_flag = 1;

  nelements = 0;
  elements = NULL;
  nparams = maxparam = 0;
  params = NULL;
  elem2param = NULL;

  nmax = 0;

  pgsize = oneatom = 0;

  calc_status = CALC_FINISHED;
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairTersoff_k2::~PairTersoff_k2()
{
  if (elements)
    for (int i = 0; i < nelements; i++) delete [] elements[i];
  delete [] elements;
  memory->destroy(params);
  memory->destroy(elem2param);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    delete [] map;
  }
}

/* ---------------------------------------------------------------------- */

void PairTersoff_k2::compute_pre_qeq(int eflag, int vflag)
{
  if(calc_status == PRE_QEQ_FINISHED){
    return;
  }

  double** const x = atom->x;
  double** const f = atom->f;
  tagint* const tag = atom->tag;
  int* const type = atom->type;

  const int inum = list->inum;
  int* const ilist = list->ilist;
  int* const numneigh = list->numneigh;
  int** const firstneigh = list->firstneigh;

  double* const special_coul = force->special_coul;

  if (atom->nmax > nmax) {
    nmax = atom->nmax;
    pre_qeq_i.resize(nmax, precalc_qeq_i());

    H.firstnbr.resize(atom->nmax, 0);
   	H.numnbrs.resize(atom->nmax, 0);

    Sht.sht_top.resize(nmax, 0);
    Sht.sht_num.resize(nmax, 0);
  }

  H.jlist.clear();
  H.val.clear();
  H.fdivqiqj.clear();

  Sht.ipage.clear();
  Sht.ijpage.clear();
  // loop over full neighbor list of my atoms

	int nshort_count = 0;
  for (int ii = 0; ii < inum; ii++) {
    int i = ilist[ii];
    i &= NEIGHMASK;
    const tagint itag = tag[i];
    const int itype = map[type[i]];
    const double xtmp = x[i][0];
    const double ytmp = x[i][1];
    const double ztmp = x[i][2];

    const int iparam_i = elem2param[itype][itype][itype];
    pre_qeq_i[i].chi = params[iparam_i].el_chi;
    pre_qeq_i[i].bigj_ii = params[iparam_i].el_eta;
    const double e_shift = erfc(params[iparam_i].w_alf*params[iparam_i].coulcut)/params[iparam_i].coulcut;
    const double e_self_wolf_divqiqi = -(e_shift/2.0 + params[iparam_i].w_alf/MY_PIS) *force->qqrd2e;
    pre_qeq_i[i].bigj_ii += 2.0*e_self_wolf_divqiqi;

    pre_qeq_i[i].bigj_ii += params[iparam_i].gamma_ecoeff*force->qqrd2e;

    int* const jlist = firstneigh[i];
    const int jnum = numneigh[i];

    H.firstnbr[i] = H.jlist.size();
    Sht.sht_top[i] = Sht.ipage.size();

    // two-body interactions

		int nshortj = 0;
    for (int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj];
      j &= NEIGHMASK;
      const tagint jtag = tag[j];

      const double delx = xtmp - x[j][0];
      const double dely = ytmp - x[j][1];
      const double delz = ztmp - x[j][2];
      const double rsq = delx*delx + dely*dely + delz*delz;

      const int jtype = map[type[j]];
      const int iparam_ij = elem2param[itype][jtype][jtype];

      if(rsq > params[iparam_ij].coulcutsq) continue;

      bool long_flag = true;
      if (itag > jtag) {
        if ((itag+jtag) % 2 == 0) long_flag = false;
      } else if (itag < jtag) {
        if ((itag+jtag) % 2 == 1) long_flag = false;
      } else {
        if (x[j][2] < x[i][2]) long_flag = false;
        if (x[j][2] == ztmp && x[j][1] < ytmp) long_flag = false;
        if (x[j][2] == ztmp && x[j][1] == ytmp && x[j][0] < xtmp) long_flag = false;
      }

      if(long_flag){

       double fvionij_divqij, vionij_divqij;
       coulomb_interaction(&params[iparam_ij],rsq,special_coul[sbmask(j)],fvionij_divqij,1,vionij_divqij);
       H.jlist.push_back(j);
       H.val.push_back(vionij_divqij);
       H.fdivqiqj.push_back(fvionij_divqij);
      }

      if (rsq > params[iparam_ij].cutsq) continue;

      Sht.ipage.push_back(j);
      nshortj++;
      nshort_count++;
    }
    H.numnbrs[i] = H.jlist.size() - H.firstnbr[i];
    Sht.sht_num[i] = nshortj;
	}

  Sht.ijpage.resize(nshort_count, precalc_qeq_ij());

  for (int ii = 0; ii < inum; ii++) {
    int i = ilist[ii];
    const tagint itag = tag[i];
    const int itype = map[type[i]];
    const double xtmp = x[i][0];
    const double ytmp = x[i][1];
    const double ztmp = x[i][2];
    const int sht_first = Sht.sht_top[i];
    const int sht_jnum = Sht.sht_num[i];

    // three-body interactions
    // skip immediately if I-J is not within cutoff

    for (int jj = 0; jj < sht_jnum; jj++) {
      int j = Sht.ipage[sht_first+jj];
      j &= NEIGHMASK;
      const int jtype = map[type[j]];
      const int iparam_ij = elem2param[itype][jtype][jtype];

      double delr1[3];
      delr1[0] = x[j][0] - xtmp;
      delr1[1] = x[j][1] - ytmp;
      delr1[2] = x[j][2] - ztmp;
      const double rsq1 = delr1[0]*delr1[0] + delr1[1]*delr1[1] + delr1[2]*delr1[2];
      if (rsq1 > params[iparam_ij].cutsq_ij) continue;

      // accumulate bondorder zeta for each i-j interaction via loop over k

      double zeta_ij1 = 0.0;
      double zeta_ij2 = 0.0;

      for (int kk = 0; kk < sht_jnum; kk++) {
        if (jj == kk) continue;
        int k = Sht.ipage[sht_first+kk];
        k &= NEIGHMASK;
        const int ktype = map[type[k]];
        const int iparam_ijk = elem2param[itype][jtype][ktype];

        double delr2[3];
        delr2[0] = x[k][0] - xtmp;
        delr2[1] = x[k][1] - ytmp;
        delr2[2] = x[k][2] - ztmp;
        const double rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];
        if (rsq2 > params[iparam_ijk].cutsq_ik) continue;

        double zeta1, zeta2;
        zeta(&params[iparam_ijk],rsq1,rsq2,delr1,delr2, zeta1, zeta2);
        zeta_ij1 += zeta1;
        zeta_ij2 += zeta2;
      }

      double prefactor1, prefactor2, fpair_divfq, evdwl_divfq;
      force_zeta(&params[iparam_ij],rsq1,zeta_ij1,zeta_ij2,fpair_divfq,prefactor1,prefactor2,1,evdwl_divfq);
      Sht.ijpage[sht_first+jj].prefactor1 = prefactor1;
      Sht.ijpage[sht_first+jj].prefactor2 = prefactor2;
      Sht.ijpage[sht_first+jj].fpair_divfq = fpair_divfq;
      Sht.ijpage[sht_first+jj].evdwl_divfq = evdwl_divfq;
      pre_qeq_i[i].chi += ters_fq_d(&params[iparam_ij],0)*evdwl_divfq;
      pre_qeq_i[i].bigj_ii += 2.0*0.5*ters_fq_2d(&params[iparam_ij],0)*evdwl_divfq;
    }
  }

  calc_status = PRE_QEQ_FINISHED;
}

void PairTersoff_k2::compute(int eflag, int vflag)
{
  compute_pre_qeq(eflag, vflag);

  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = vflag_atom = 0;

  double fpair, evdwl, ecoul;
  double **x = atom->x;
  double **f = atom->f;
  double* q = atom->q;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  tagint *tag = atom->tag;
  int inum = list->inum;
  int* ilist = list->ilist;

  double *special_coul = force->special_coul;

  // loop over full neighbor list of my atoms

  for (int ii = 0; ii < inum; ii++) {
    int i = ilist[ii];
    tagint itag = tag[i];
    int itype = map[type[i]];
    double xtmp = x[i][0];
    double ytmp = x[i][1];
    double ztmp = x[i][2];

    int iparam_i = elem2param[itype][itype][itype];
    double e_shift = erfc(params[iparam_i].w_alf*params[iparam_i].coulcut)/params[iparam_i].coulcut;
    double e_self_wolf = -(e_shift/2.0 + params[iparam_i].w_alf/MY_PIS) * q[i]*q[i]*force->qqrd2e;
    double e_self_chi = (params[iparam_i].el_chi*q[i]+0.5*params[iparam_i].el_eta*q[i]*q[i]);
    double e_self = e_self_wolf+e_self_chi+params[iparam_i].u_shift;
    if (evflag) ev_tally(i,i,nlocal,0,0.0,e_self,0.0,0.0,0.0,0.0);

    // two-body interactions, skip half of them
    for(int jj = H.firstnbr[i]; jj < H.firstnbr[i]+H.numnbrs[i]; jj++){
      int j = H.jlist[jj];
      double delx = xtmp - x[j][0];
      double dely = ytmp - x[j][1];
      double delz = ztmp - x[j][2];

      double fvionij = H.fdivqiqj[jj]*q[i]*q[j];
      double vionij = H.val[jj]*q[i]*q[j];

      f[i][0] += delx*fvionij;
      f[i][1] += dely*fvionij;
      f[i][2] += delz*fvionij;
      f[j][0] -= delx*fvionij;
      f[j][1] -= dely*fvionij;
      f[j][2] -= delz*fvionij;

      if(evflag) ev_tally(i,j,nlocal,newton_pair,0.0,vionij,fvionij,delx,dely,delz);
    }

    int sht_first = Sht.sht_top[i];
    int sht_jnum = Sht.sht_num[i];

    for (int jj = 0; jj < sht_jnum; jj++) {
      int j = Sht.ipage[sht_first+jj];
      j &= NEIGHMASK;
      tagint jtag = tag[j];

      if (itag > jtag) {
        if ((itag+jtag) % 2 == 0) continue;
      } else if (itag < jtag) {
        if ((itag+jtag) % 2 == 1) continue;
      } else {
        if (x[j][2] < x[i][2]) continue;
        if (x[j][2] == ztmp && x[j][1] < ytmp) continue;
        if (x[j][2] == ztmp && x[j][1] == ytmp && x[j][0] < xtmp) continue;
      }

      int jtype = map[type[j]];

      double delx = xtmp - x[j][0];
      double dely = ytmp - x[j][1];
      double delz = ztmp - x[j][2];
      double rsq = delx*delx + dely*dely + delz*delz;

      int iparam_ij = elem2param[itype][jtype][jtype];

      if (rsq > params[iparam_ij].cutsq_ij) continue;

      repulsive(&params[iparam_ij],rsq,fpair,eflag,evdwl);

      f[i][0] += delx*fpair;
      f[i][1] += dely*fpair;
      f[i][2] += delz*fpair;
      f[j][0] -= delx*fpair;
      f[j][1] -= dely*fpair;
      f[j][2] -= delz*fpair;

      if (evflag) ev_tally(i,j,nlocal,newton_pair,
                           evdwl,0.0,fpair,delx,dely,delz);
    }
    // three-body interactions
    // skip immediately if I-J is not within cutoff

    for (int jj = 0; jj < sht_jnum; jj++) {
      int j = Sht.ipage[sht_first+jj];
      j &= NEIGHMASK;
      int jtype = map[type[j]];
      int iparam_ij = elem2param[itype][jtype][jtype];

      double delr1[3];
      delr1[0] = x[j][0] - xtmp;
      delr1[1] = x[j][1] - ytmp;
      delr1[2] = x[j][2] - ztmp;
      double rsq1 = delr1[0]*delr1[0] + delr1[1]*delr1[1] + delr1[2]*delr1[2];
      if (rsq1 > params[iparam_ij].cutsq_ij) continue;

      // pairwise force due to zeta

      double fq = ters_fq(&params[iparam_ij], q[i]);
      double fq_0 = ters_fq(&params[iparam_ij], 0);
      fpair = Sht.ijpage[sht_first+jj].fpair_divfq*fq;
      evdwl = Sht.ijpage[sht_first+jj].evdwl_divfq*fq_0;
      ecoul = Sht.ijpage[sht_first+jj].evdwl_divfq*(fq-fq_0);

      f[i][0] += delr1[0]*fpair;
      f[i][1] += delr1[1]*fpair;
      f[i][2] += delr1[2]*fpair;
      f[j][0] -= delr1[0]*fpair;
      f[j][1] -= delr1[1]*fpair;
      f[j][2] -= delr1[2]*fpair;

      if (evflag) ev_tally(i,j,nlocal,newton_pair,
                           evdwl,ecoul,-fpair,-delr1[0],-delr1[1],-delr1[2]);

      // attractive term via loop over k

      for (int kk = 0; kk < sht_jnum; kk++) {
        if (jj == kk) continue;
        int k = Sht.ipage[sht_first+kk];
        k &= NEIGHMASK;
        int ktype = map[type[k]];
        int iparam_ijk = elem2param[itype][jtype][ktype];

        double delr2[3];
        delr2[0] = x[k][0] - xtmp;
        delr2[1] = x[k][1] - ytmp;
        delr2[2] = x[k][2] - ztmp;
        double rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];
        if (rsq2 > params[iparam_ijk].cutsq_ik) continue;

        double fi[3],fj[3],fk[3];
        double prefactor1 = Sht.ijpage[sht_first+jj].prefactor1;
        prefactor1 *= ters_fq(&params[iparam_ijk], q[i]);
        attractive(&params[iparam_ijk],params[iparam_ijk].powermint,params[iparam_ijk].alpha,
            params[iparam_ijk].h,params[iparam_ijk].c,params[iparam_ijk].d,params[iparam_ijk].bigre_ij-params[iparam_ijk].bigre_ik,prefactor1,
                   rsq1,rsq2,delr1,delr2,fi,fj,fk);

        f[i][0] += fi[0];
        f[i][1] += fi[1];
        f[i][2] += fi[2];
        f[j][0] += fj[0];
        f[j][1] += fj[1];
        f[j][2] += fj[2];
        f[k][0] += fk[0];
        f[k][1] += fk[1];
        f[k][2] += fk[2];

        if (vflag_atom) v_tally3(i,j,k,fj,fk,delr1,delr2);

        double prefactor2 = Sht.ijpage[sht_first+jj].prefactor2;
        prefactor2 *= ters_fq(&params[iparam_ijk], q[i]);
        attractive(&params[iparam_ijk],params[iparam_ijk].powermint2,params[iparam_ijk].alpha2,
            params[iparam_ijk].h2,params[iparam_ijk].c2,params[iparam_ijk].d2,params[iparam_ijk].bigre_ij2-params[iparam_ijk].bigre_ik2,prefactor2,
                   rsq1,rsq2,delr1,delr2,fi,fj,fk);

        f[i][0] += fi[0];
        f[i][1] += fi[1];
        f[i][2] += fi[2];
        f[j][0] += fj[0];
        f[j][1] += fj[1];
        f[j][2] += fj[2];
        f[k][0] += fk[0];
        f[k][1] += fk[1];
        f[k][2] += fk[2];

        if (vflag_atom) v_tally3(i,j,k,fj,fk,delr1,delr2);
      }
    }
  }
  if (vflag_fdotr) virial_fdotr_compute();

  calc_status = CALC_FINISHED;
}

/* ---------------------------------------------------------------------- */

void PairTersoff_k2::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  map = new int[n+1];
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairTersoff_k2::settings(int narg, char **arg)
{
  if (narg != 0) error->all(FLERR,"Illegal pair_style command");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairTersoff_k2::coeff(int narg, char **arg)
{
  int i,j,n;

  if (!allocated) allocate();

  if (narg != 3 + atom->ntypes)
    error->all(FLERR,"Incorrect args for pair coefficients");

  // insure I,J args are * *

  if (strcmp(arg[0],"*") != 0 || strcmp(arg[1],"*") != 0)
    error->all(FLERR,"Incorrect args for pair coefficients");

  // read args that map atom types to elements in potential file
  // map[i] = which element the Ith atom type is, -1 if NULL
  // nelements = # of unique elements
  // elements = list of element names

  if (elements) {
    for (i = 0; i < nelements; i++) delete [] elements[i];
    delete [] elements;
  }
  elements = new char*[atom->ntypes];
  for (i = 0; i < atom->ntypes; i++){
    elements[i] = NULL;
  }

  nelements = 0;
  for (i = 3; i < narg; i++) {
    if (strcmp(arg[i],"NULL") == 0) {
      map[i-2] = -1;
      continue;
    }
    for (j = 0; j < nelements; j++)
      if (strcmp(arg[i],elements[j]) == 0) break;
    map[i-2] = j;
    if (j == nelements) {
      n = strlen(arg[i]) + 1;
      elements[j] = new char[n];
      strcpy(elements[j],arg[i]);
      nelements++;
    }
  }

  // read potential file and initialize potential parameters

  read_file(arg[2]);

  // clear setflag since coeff() called once with I,J = * *

  n = atom->ntypes;
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  // set setflag i,j for type pairs where both are mapped to elements

  int count = 0;
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      if (map[i] >= 0 && map[j] >= 0) {
        setflag[i][j] = 1;
        count++;
      }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairTersoff_k2::init_style()
{
  if (atom->tag_enable == 0)
    error->all(FLERR,"Pair style Tersoff requires atom IDs");
  if (force->newton_pair == 0)
    error->all(FLERR,"Pair style Tersoff requires newton pair on");
  if (!atom->q_flag)
    error->all(FLERR,"Pair style Tersoff/k requires atom attribute q");

  // need a full neighbor list

  int irequest = neighbor->request(this);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;

  // local Comb neighbor list
  // create pages if first time or if neighbor pgsize/oneatom has changed

  int create = 0;
  if (pgsize != neighbor->pgsize) create = 1;
  if (oneatom != neighbor->oneatom) create = 1;

  if (create) {
    pgsize = neighbor->pgsize;
    oneatom = neighbor->oneatom;
  }
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairTersoff_k2::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  return cutmax;
}

/* ---------------------------------------------------------------------- */

void PairTersoff_k2::read_file(char *file)
{
  const int params_per_line_max = 40;
  char *words[params_per_line_max+1];

  memory->sfree(params);
  params = NULL;
  nparams = maxparam = 0;

  // open file on proc 0

  FILE *fp;
  if (comm->me == 0) {
    fp = force->open_potential(file);
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open Tersoff potential file %s",file);
      error->one(FLERR,str);
    }
  }

  // read each line out of file, skipping blank lines or leading '#'
  // store line of params if all 3 element tags are in element list

  int n,nwords,ielement,jelement,kelement;
  char line[MAXLINE],*ptr;
  int eof = 0;
  enum param_type_enum{
    SINGLE,DOUBLE,TRIPLE
  } param_type;
  std::vector<param_single_struct> params_single;
  std::vector<param_double_struct> params_double;

  while (1) {
    if (comm->me == 0) {
      ptr = fgets(line,MAXLINE,fp);
      if (ptr == NULL) {
        eof = 1;
        fclose(fp);
      } else n = strlen(line) + 1;
    }
    MPI_Bcast(&eof,1,MPI_INT,0,world);
    if (eof) break;
    MPI_Bcast(&n,1,MPI_INT,0,world);
    MPI_Bcast(line,n,MPI_CHAR,0,world);

    // strip comment, skip line if blank

    if ((ptr = strchr(line,'#'))) *ptr = '\0';
    nwords = atom->count_words(line);
    if (nwords == 0) continue;

    // concatenate additional lines until have params_per_line words

    int params_per_line;
    if(strncmp(line, "single", 6) == 0){
      param_type = SINGLE;
      params_per_line = 7;
    }else if(strncmp(line, "double", 6) == 0){
      param_type = DOUBLE;
      params_per_line = 35;
    }else if(strncmp(line, "triple", 6) == 0){
      param_type = TRIPLE;
      params_per_line = 14;
    }else{
      error->all(FLERR,"Incorrect format in Tersoff/k potential file (no single/double/triple word)");
    }

    while (nwords < params_per_line) {
      n = strlen(line);
      if (comm->me == 0) {
        ptr = fgets(&line[n],MAXLINE-n,fp);
        if (ptr == NULL) {
          eof = 1;
          fclose(fp);
        } else n = strlen(line) + 1;
      }
      MPI_Bcast(&eof,1,MPI_INT,0,world);
      if (eof) break;
      MPI_Bcast(&n,1,MPI_INT,0,world);
      MPI_Bcast(line,n,MPI_CHAR,0,world);
      if ((ptr = strchr(line,'#'))) *ptr = '\0';
      nwords = atom->count_words(line);
    }

    if (nwords != params_per_line){
      error->all(FLERR,"Incorrect format in Tersoff/k potential file");
    }

    // words = ptrs to all words in line

    nwords = 0;
    words[nwords++] = strtok(line," \t\n\r\f");
    while ((words[nwords++] = strtok(NULL," \t\n\r\f"))) continue;

    switch(param_type){
      case SINGLE:
        {
          for(ielement = 0; ielement < nelements; ielement++)
            if(strcmp(words[1],elements[ielement]) == 0) break;
          if(ielement == nelements) continue;
          param_single_struct tmp_single;
          tmp_single.ielement = ielement;
          tmp_single.u_shift = atof(words[2]);
          tmp_single.el_chi = atof(words[3]);
          tmp_single.el_eta = atof(words[4]);
          tmp_single.ni_newtral = atof(words[5]);
          tmp_single.ni0 = atof(words[6]);
          params_single.push_back(tmp_single);
        }
        break;
      case DOUBLE:
        {
          for (ielement = 0; ielement < nelements; ielement++)
            if (strcmp(words[1],elements[ielement]) == 0) break;
          if (ielement == nelements) continue;
          for (jelement = 0; jelement < nelements; jelement++)
            if (strcmp(words[2],elements[jelement]) == 0) break;
          if (jelement == nelements) continue;
          param_double_struct tmp_double;
          tmp_double.ielement = ielement;
          tmp_double.jelement = jelement;
          tmp_double.biga1 = atof(words[3]);
          tmp_double.biga2 = atof(words[4]);
          tmp_double.biga3 = atof(words[5]);
          tmp_double.bigb1 = atof(words[6]);
          tmp_double.bigb2 = atof(words[7]);
          tmp_double.bigb3 = atof(words[8]);
          tmp_double.lama1 = atof(words[9]);
          tmp_double.lama2 = atof(words[10]);
          tmp_double.lama3 = atof(words[11]);
          tmp_double.lamb1 = atof(words[12]);
          tmp_double.lamb2 = atof(words[13]);
          tmp_double.lamb3 = atof(words[14]);
          tmp_double.powern = atof(words[15]);
          tmp_double.powern2 = atof(words[16]);
          tmp_double.powern_del = atof(words[17]);
          tmp_double.powern_rev = atof(words[18]);
          tmp_double.powern2_rev = atof(words[19]);
          tmp_double.powern_rev_del = atof(words[20]);
          tmp_double.powerp = atof(words[21]);
          tmp_double.powerp_rev = atof(words[22]);
          tmp_double.zeta_shift1 = atof(words[23]);
          tmp_double.zeta_shift2 = atof(words[24]);
          tmp_double.zeta_shift1_rev = atof(words[25]);
          tmp_double.zeta_shift2_rev = atof(words[26]);
          tmp_double.Q = atof(words[27]);
          tmp_double.bigre = atof(words[28]);
          tmp_double.bigre2 = atof(words[29]);
          tmp_double.bigr = atof(words[30]);
          tmp_double.bigd = atof(words[31]);
          tmp_double.w_alpha = atof(words[32]);
          tmp_double.coulcut = atof(words[33]);
          tmp_double.gamma = atof(words[34]);
          params_double.push_back(tmp_double);
        }
        break;
      case TRIPLE:
        {

          // ielement,jelement,kelement = 1st args
          // if all 3 args are in element list, then parse this line
          // else skip to next line

          for (ielement = 0; ielement < nelements; ielement++)
            if (strcmp(words[1],elements[ielement]) == 0) break;
          if (ielement == nelements) continue;
          for (jelement = 0; jelement < nelements; jelement++)
            if (strcmp(words[2],elements[jelement]) == 0) break;
          if (jelement == nelements) continue;
          for (kelement = 0; kelement < nelements; kelement++)
            if (strcmp(words[3],elements[kelement]) == 0) break;
          if (kelement == nelements) continue;

          // load up parameter settings and error check their values

          if (nparams == maxparam) {
            maxparam += DELTA;
            params = (Param *) memory->srealloc(params,maxparam*sizeof(Param),
                "pair:params");
          }
          params[nparams].ielement = ielement;
          params[nparams].jelement = jelement;
          params[nparams].kelement = kelement;
          params[nparams].powerm = atof(words[4]);
          params[nparams].alpha = atof(words[5]);
          params[nparams].c = atof(words[6]);
          params[nparams].d = atof(words[7]);
          params[nparams].h = atof(words[8]);
          params[nparams].powerm2 = atof(words[9]);
          params[nparams].alpha2 = atof(words[10]);
          params[nparams].c2 = atof(words[11]);
          params[nparams].d2 = atof(words[12]);
          params[nparams].h2 = atof(words[13]);

          // currently only allow m exponent of 1 or 3
          params[nparams].powermint = int(params[nparams].powerm);
          params[nparams].powermint2 = int(params[nparams].powerm2);
          /*
             if (params[nparams].c < 0.0 || params[nparams].d < 0.0 ||
             params[nparams].powern < 0.0 || params[nparams].beta < 0.0 ||
             params[nparams].lam2 < 0.0 || params[nparams].bigb < 0.0 ||
             params[nparams].bigr < 0.0 ||params[nparams].bigd < 0.0 ||
             params[nparams].bigd > params[nparams].bigr ||
             params[nparams].lam1 < 0.0 || params[nparams].biga < 0.0 ||
             params[nparams].powerm - params[nparams].powermint != 0.0 ||
             (params[nparams].powermint != 3 && params[nparams].powermint != 1) ||
             params[nparams].gamma < 0.0)
             error->all(FLERR,"Illegal Tersoff parameter");
           */
          nparams++;
        }
        break;
    }
  }

  // set elem2param for all element triplet combinations
  // must be a single exact match to lines read from file
  // do not allow for ACB in place of ABC

  memory->destroy(elem2param);
  memory->create(elem2param,nelements,nelements,nelements,"pair:elem2param");

  std::vector<int> elem2param_single;
  std::vector<std::vector<int> > elem2param_double;

  for(int i = 0; i < nelements; i++){
    elem2param_single.push_back(-1);
    elem2param_double.push_back(std::vector<int>());
    for(int j = 0; j < nelements; j++){
      elem2param_double[i].push_back(-1);
      for(int k = 0; k < nelements; k++){
        elem2param[i][j][k] = -1;
      }
    }
  }

  for(unsigned int i = 0; i < params_single.size(); i++){
    int pi = params_single[i].ielement;
    if(elem2param_single[pi] != -1){
      error->all(FLERR,"Potential file has duplicate entry(single)");
    }
    elem2param_single[pi] = i;
  }

  for(unsigned int i = 0; i < params_double.size(); i++){
    int pi = params_double[i].ielement;
    int pj = params_double[i].jelement;
    if(elem2param_double[pi][pj] != -1){
      error->all(FLERR,"Potential file has duplicate entry(double)");
    }
    elem2param_double[pi][pj] = i;
    elem2param_double[pj][pi] = i;
  }

  for(int i = 0; i < nparams; i++){
    int pi = params[i].ielement;
    int pj = params[i].jelement;
    int pk = params[i].kelement;
    if(elem2param[pi][pj][pk] != -1){
      error->all(FLERR,"Potential file has duplicate entry(triple)");
    }
    elem2param[pi][pj][pk] = i;
  }

  for(int i = 0; i < nelements; i++){
    if(elem2param_single[i] == -1){
      error->all(FLERR,"Potential file is missing an entry(single)");
    }
    int iparam_i = elem2param_single[i];
    for(int j = 0; j < nelements; j++){
      if(elem2param_double[i][j] == -1){
        error->all(FLERR,"Potential file is missing an entry(double)");
      }
      int iparam_ij = elem2param_double[i][j];
      for(int k = 0; k < nelements; k++){
        int iparam_ik = elem2param_double[i][k];
        if(elem2param[i][j][k] == -1){
          error->all(FLERR,"Potential file is missing an entry(triple)");
        }
        int iparam_ijk = elem2param[i][j][k];
        params[iparam_ijk].u_shift = params_single[iparam_i].u_shift;
        params[iparam_ijk].el_chi = params_single[iparam_i].el_chi;
        params[iparam_ijk].el_eta = params_single[iparam_i].el_eta;
        params[iparam_ijk].ni_newtral = params_single[iparam_i].ni_newtral;
        params[iparam_ijk].ni0 = params_single[iparam_i].ni0;
        params[iparam_ijk].biga1 = params_double[iparam_ij].biga1;
        params[iparam_ijk].biga2 = params_double[iparam_ij].biga2;
        params[iparam_ijk].biga3 = params_double[iparam_ij].biga3;
        params[iparam_ijk].bigb1 = params_double[iparam_ij].bigb1;
        params[iparam_ijk].bigb2 = params_double[iparam_ij].bigb2;
        params[iparam_ijk].bigb3 = params_double[iparam_ij].bigb3;
        params[iparam_ijk].lama1 = params_double[iparam_ij].lama1;
        params[iparam_ijk].lama2 = params_double[iparam_ij].lama2;
        params[iparam_ijk].lama3 = params_double[iparam_ij].lama3;
        params[iparam_ijk].lamb1 = params_double[iparam_ij].lamb1;
        params[iparam_ijk].lamb2 = params_double[iparam_ij].lamb2;
        params[iparam_ijk].lamb3 = params_double[iparam_ij].lamb3;
        if(params[iparam_ijk].ielement == params_double[iparam_ij].ielement && params[iparam_ijk].jelement == params_double[iparam_ij].jelement){
          params[iparam_ijk].powern = params_double[iparam_ij].powern;
          params[iparam_ijk].powern2 = params_double[iparam_ij].powern2;
          params[iparam_ijk].powern_del = params_double[iparam_ij].powern_del;
          params[iparam_ijk].powerp = params_double[iparam_ij].powerp;
          params[iparam_ijk].zeta_shift1 = params_double[iparam_ij].zeta_shift1;
          params[iparam_ijk].zeta_shift2 = params_double[iparam_ij].zeta_shift2;
        }else{
          params[iparam_ijk].powern = params_double[iparam_ij].powern_rev;
          params[iparam_ijk].powern2 = params_double[iparam_ij].powern2_rev;
          params[iparam_ijk].powern_del = params_double[iparam_ij].powern_rev_del;
          params[iparam_ijk].powerp = params_double[iparam_ij].powerp_rev;
          params[iparam_ijk].zeta_shift1 = params_double[iparam_ij].zeta_shift1_rev;
          params[iparam_ijk].zeta_shift2 = params_double[iparam_ij].zeta_shift2_rev;
        }
        params[iparam_ijk].zeta_shift_total = pow(pow(params[iparam_ijk].zeta_shift1, -1.0*params[iparam_ijk].powerp)+pow(params[iparam_ijk].zeta_shift1, -1.0*params[iparam_ijk].powerp), -1.0/(2.0*params[iparam_ijk].powern_del*params[iparam_ijk].powerp));
        params[iparam_ijk].Q = params_double[iparam_ij].Q;
        params[iparam_ijk].bigre_ij = params_double[iparam_ij].bigre;
        params[iparam_ijk].bigre_ik = params_double[iparam_ik].bigre;
        params[iparam_ijk].bigre_ij2 = params_double[iparam_ij].bigre2;
        params[iparam_ijk].bigre_ik2 = params_double[iparam_ik].bigre2;
        params[iparam_ijk].bigr_ij = params_double[iparam_ij].bigr;
        params[iparam_ijk].bigr_ik = params_double[iparam_ik].bigr;
        params[iparam_ijk].bigd_ij = params_double[iparam_ij].bigd;
        params[iparam_ijk].bigd_ik = params_double[iparam_ik].bigd;
        params[iparam_ijk].w_alf = params_double[iparam_ij].w_alpha;
        params[iparam_ijk].coulcut = params_double[iparam_ij].coulcut;
        params[iparam_ijk].gamma = params_double[iparam_ij].gamma;
      }
    }
  }


  // compute parameter values derived from inputs

  for (int m = 0; m < nparams; m++) {
    params[m].cutoff_corr_ij = log(0.9)/(1.0/(params[m].bigr_ij+params[m].bigd_ij)-1.0/(2.0*params[m].bigd_ij));
    params[m].cutoff_r0_ij = exp(-params[m].cutoff_corr_ij/(params[m].bigr_ij+params[m].bigd_ij));
    params[m].cutoff_corr_ik = log(0.9)/(1.0/(params[m].bigr_ik+params[m].bigd_ik)-1.0/(2.0*params[m].bigd_ik));
    params[m].cutoff_r0_ik = exp(-params[m].cutoff_corr_ij/(params[m].bigr_ik+params[m].bigd_ik));

		params[m].cutsq_ij = (params[m].bigr_ij + params[m].bigd_ij)*(params[m].bigr_ij + params[m].bigd_ij);
		params[m].cutsq_ik = (params[m].bigr_ik + params[m].bigd_ik)*(params[m].bigr_ik + params[m].bigd_ik);
    params[m].cut = fmax(params[m].bigr_ij + params[m].bigd_ij, params[m].bigr_ik+ params[m].bigd_ik);
    params[m].cutsq = params[m].cut*params[m].cut;
    params[m].coulcutsq = params[m].coulcut*params[m].coulcut;

    params[m].ca1 = pow(2.0*params[m].powern_del*1.0e-16,-1.0/params[m].powern);
    params[m].ca1_2 = pow(2.0*params[m].powern_del*1.0e-16,-1.0/params[m].powern2);
    params[m].ca4 = 1.0/params[m].ca1;

    const double cut_coul = params[m].coulcut;
    const double alf = params[m].w_alf;
    params[m].e_shift = erfc(alf*cut_coul)/cut_coul;
    params[m].f_shift = -(params[m].e_shift+ 2.0*alf/MY_PIS * exp(-alf*alf*cut_coul*cut_coul)) / cut_coul;
    const double gamma = params[m].gamma;
    const double r_effc = pow(cut_coul*cut_coul*cut_coul+(1.0/gamma/gamma/gamma), 1.0/3.0);
    params[m].gamma_ecoeff = 1.0/r_effc - 1.0/cut_coul;
    params[m].gamma_fcoeff = cut_coul*cut_coul/r_effc/r_effc/r_effc/r_effc-1.0/cut_coul/cut_coul;
  }

  // set cutmax to max of all params

  cutmax = 0.0;
  cutshort = 0.0;
  for (int m = 0; m < nparams; m++){
    if (params[m].cut > cutmax) cutmax = params[m].cut;
    if (params[m].coulcut > cutmax) cutmax = params[m].coulcut;
    if (params[m].cut > cutshort) cutshort = params[m].cut;
  }
}

/* ---------------------------------------------------------------------- */

void PairTersoff_k2::coulomb_interaction(const Param* const param, const double rsq, const double factor_coul, double &fvionij_divq, const int eflag, double& eng_divq){
  const double alf = param->w_alf;
  const double e_shift = param->e_shift;//erfc(alf*cut_coul)/cut_coul;
  const double f_shift = param->f_shift;//-(e_shift+ 2.0*alf/MY_PIS * exp(-alf*alf*cut_coul*cut_coul)) / cut_coul;
  const double r = sqrt(rsq);
  const double prefactor = force->qqrd2e/r;

  const double gamma = param->gamma;
  const double r_eff = pow(r*r*r+1.0/gamma/gamma/gamma, 1.0/3.0);
  const double fvslater_divq = (rsq/r_eff/r_eff/r_eff/r_eff-1.0/rsq-param->gamma_fcoeff);

  const double erfcc = erfc(alf*r);
  const double erfcd = exp(-alf*alf*rsq);
  const double dvdrr = erfcc/rsq + 2.0*alf/MY_PIS * erfcd/r + f_shift;
  fvionij_divq = (dvdrr+fvslater_divq)*prefactor;
  if (factor_coul < 1.0) fvionij_divq -= (1.0-factor_coul)*prefactor/rsq;

  if (eflag) {
    const double v_slater = (r/r_eff-1.0-r*param->gamma_ecoeff)*prefactor;
    const double v_sh = (erfcc - e_shift*r) * prefactor;
    eng_divq = v_slater + v_sh + (param->gamma_fcoeff - f_shift)*(r-param->coulcut)*r*prefactor;
    if (factor_coul < 1.0) eng_divq -= (1.0-factor_coul)*prefactor;
  }
}


/* ---------------------------------------------------------------------- */

void PairTersoff_k2::repulsive(Param* const param, const double rsq, double &fforce,
                            const int eflag, double &eng)
{
  const double tmp_fQ = k_fQ(param, rsq);
  const double tmp_fQ_d = k_fQ_d(param, rsq);

  const double r = sqrt(rsq);
  const double tmp_fc = ters_fc(0, r,param);
  const double tmp_fc_d = ters_fc_d(0, r,param);
  const double tmp_exp1 = param->biga1 * exp(-param->lama1 * r);
  const double tmp_exp2 = param->biga2 * exp(-param->lama2 * r);
  const double tmp_exp3 = param->biga3 * exp(-param->lama3 * r);
  const double tmp_exp = tmp_exp1+tmp_exp2+tmp_exp3;
  const double tmp_exp_d = param->lama1*tmp_exp1+param->lama2*tmp_exp2+param->lama3*tmp_exp3;
  fforce = -1.0*(tmp_fc_d*tmp_fQ*tmp_exp + tmp_fc*tmp_fQ_d*tmp_exp - tmp_fc*tmp_fQ*(tmp_exp_d)) / r;
  if (eflag) eng = tmp_fc * tmp_fQ * tmp_exp;
}

/* ---------------------------------------------------------------------- */

void PairTersoff_k2::zeta(Param* const param, const double rsqij, const double rsqik,
                         double *delrij, double *delrik, double& zeta1, double& zeta2)
{
  double costheta,arg1,arg2,ex_delr1,ex_delr2;

  const double rij = sqrt(rsqij);
  const double rik = sqrt(rsqik);
  costheta = (delrij[0]*delrik[0] + delrij[1]*delrik[1] +
              delrij[2]*delrik[2]) / (rij*rik);

	double rdiff1 = (rij-param->bigre_ij)-(rik-param->bigre_ik);
	double rdiff2 = (rij-param->bigre_ij2)-(rik-param->bigre_ik2);
  if (param->powermint == 3){
    arg1 = param->alpha*pow(rdiff1,3.0);
  }else{
    arg1 = param->alpha * rdiff1;
  }
  if (param->powermint2 == 3){
    arg2 = param->alpha2*pow(rdiff2,3.0);
  }else{
    arg2 = param->alpha2 * rdiff2;
  }

  if (arg1 > 69.0776) ex_delr1 = 1.e30;
  else if (arg1 < -69.0776) ex_delr1 = 0.0;
  else ex_delr1 = exp(arg1);

  if (arg2 > 69.0776) ex_delr2 = 1.e30;
  else if (arg2 < -69.0776) ex_delr2 = 0.0;
  else ex_delr2 = exp(arg2);

  double ters_gijk1, ters_gijk2;

  zeta1 = ters_fc(1, rik,param) * ters_gijk(costheta,param->h, param->c, param->d) * ex_delr1;
  zeta2 = ters_fc(1, rik,param) * ters_gijk(costheta,param->h2, param->c2, param->d2) * ex_delr2;
}

/* ---------------------------------------------------------------------- */

void PairTersoff_k2::force_zeta(Param* const param, const double rsq, const double zeta1_ij, const double zeta2_ij, double &fforce, double &prefactor1, double &prefactor2,
                             int eflag, double &eng)
{
  const double r = sqrt(rsq);
  const double fa = ters_fa(r,param);
  const double fa_d = ters_fa_d(r,param);
  const double bij = ters_bij(zeta1_ij,zeta2_ij,param);
  double bij_d1, bij_d2;
  ters_bij_d(zeta1_ij,zeta2_ij,param,bij_d1,bij_d2);
  //const double bij_d = ters_bij_d(zeta_ij, param);
  fforce = 0.5*bij*fa_d / r;
  prefactor1 = -0.5*fa * bij_d1;
  prefactor2 = -0.5*fa * bij_d2;
  if (eflag) eng = 0.5*bij*fa;
}


/* ---------------------------------------------------------------------- */

double PairTersoff_k2::ters_fq(Param *param, const double qi)
{
  double ni_qi = param->ni_newtral-qi;
  //return 1.0;
  return ni_qi*(param->ni0 - ni_qi)/(param->ni_newtral*(param->ni0 - param->ni_newtral));
}

double PairTersoff_k2::ters_fq_d(Param *param, const double qi)
{
  double ni_qi = param->ni_newtral-qi;
  //return 0.0;
  return (2.0*ni_qi-param->ni0)/(param->ni_newtral*(param->ni0 - param->ni_newtral));
}

double PairTersoff_k2::ters_fq_2d(Param *param, const double qi)
{
  //return 0.0;
  return -2.0/(param->ni_newtral*(param->ni0 - param->ni_newtral));
}


/* ----------------------------------------------------------------------
   attractive term
   use param_ij cutoff for rij test
   use param_ijk cutoff for rik test
------------------------------------------------------------------------- */

void PairTersoff_k2::attractive(Param* const param, const double powermint, const double alpha,
                             const double h, const double c, const double d, const double bigre_diff, const double prefactor,
                             const double rsqij, const double rsqik,
                             double *delrij, double *delrik,
                             double *fi, double *fj, double *fk)
{
  double rij_hat[3],rik_hat[3];
  double rij,rijinv,rik,rikinv;

  rij = sqrt(rsqij);
  rijinv = 1.0/rij;
  vec3_scale(rijinv,delrij,rij_hat);

  rik = sqrt(rsqik);
  rikinv = 1.0/rik;
  vec3_scale(rikinv,delrik,rik_hat);

  ters_zetaterm_d(prefactor,rij_hat,rij,rik_hat,rik,fi,fj,fk,param, powermint, alpha, h, c, d, bigre_diff);
}

/* ---------------------------------------------------------------------- */

double PairTersoff_k2::ters_fc(const int ikflag, const double r, Param *param)
{
  double exp_corr = param->cutoff_corr_ij;
  double exp_r0 = param->cutoff_r0_ij;
  double ters_R = param->bigr_ij;
  double ters_D = param->bigd_ij;
  if(ikflag){
    exp_corr = param->cutoff_corr_ik;
    exp_r0 = param->cutoff_r0_ik;
    ters_R = param->bigr_ik;
    ters_D = param->bigd_ik;
  }

  double rev_r = ters_R+ters_D-r;
  if (rev_r < 0) return 0.0;
  double ans = 1.0/exp_r0*(exp(-exp_corr/(rev_r)));
  return ans;

  //if (r < ters_R-ters_D) return 1.0;
  //if (r > ters_R+ters_D) return 0.0;
	//double tmp = (r-(ters_R-ters_D))/(2.0*ters_D);
	//return (((20*tmp-70)*tmp+84)*tmp-35)*tmp*tmp*tmp*tmp+1;
  //return 0.5*(1.0 - sin(MY_PI2*(r - ters_R)/ters_D));
}

/* ---------------------------------------------------------------------- */

double PairTersoff_k2::ters_fc_d(const int ikflag, const double r, Param *param)
{
  double exp_corr = param->cutoff_corr_ij;
  double exp_r0 = param->cutoff_r0_ij;
  double ters_R = param->bigr_ij;
  double ters_D = param->bigd_ij;
  if(ikflag){
    exp_corr = param->cutoff_corr_ik;
    exp_r0 = param->cutoff_r0_ik;
    ters_R = param->bigr_ik;
    ters_D = param->bigd_ik;
  }


  if (r > ters_R+ters_D) return 0.0;
  double rev_r = ters_R+ters_D-r;
  return -1.0/exp_r0*exp_corr/(rev_r*rev_r)*(exp(-exp_corr/(rev_r)));

  //if (r < ters_R-ters_D) return 0.0;
  //if (r > ters_R+ters_D) return 0.0;
 	//double tmp = (r-(ters_R-ters_D))/(2.0*ters_D);
	//return (((140*tmp-420)*tmp+420)*tmp-140)*tmp*tmp*tmp/(2.0*ters_D);
  //return -(MY_PI4/ters_D) * cos(MY_PI2*(r - ters_R)/ters_D);
}

/* ---------------------------------------------------------------------- */

double PairTersoff_k2::ters_fa(const double r, Param *param)
{
  if (r > param->bigr_ij + param->bigd_ij) return 0.0;
  return -1.0*(param->bigb1 * exp(-param->lamb1 * r)
      + param->bigb2 * exp(-param->lamb2 * r)
      + param->bigb3 * exp(-param->lamb3 * r)
      ) * ters_fc(0, r,param);
}

/* ---------------------------------------------------------------------- */

double PairTersoff_k2::ters_fa_d(const double r, Param *param)
{
  if (r > param->bigr_ij + param->bigd_ij) return 0.0;
  double tmp_fc = ters_fc(0, r,param);
  double tmp_fc_d = ters_fc_d(0, r,param);
  return (param->bigb1 * exp(-param->lamb1 * r) *
    (param->lamb1 * tmp_fc - tmp_fc_d))
    + (param->bigb2 * exp(-param->lamb2 * r) *
    (param->lamb2 * tmp_fc - tmp_fc_d))
    + (param->bigb3 * exp(-param->lamb3 * r) *
    (param->lamb3 * tmp_fc - tmp_fc_d));
}

/* ---------------------------------------------------------------------- */

double PairTersoff_k2::ters_bij(const double zeta1, const double zeta2, Param *param)
{
  if (zeta1 > param->ca1 && zeta2 > param->ca1_2){
    if(zeta1>zeta2){
      return pow(zeta2, -param->powern2/(2.0*param->powern_del));
    }else{
      return pow(zeta1, -param->powern/(2.0*param->powern_del));
    }
    //return pow(zeta, -param->powern/(2.0*param->powern_del));
  }
  if (zeta1 < param->ca4 || zeta2 < param->ca4) return 1.0;
  const double zeta1n = pow(zeta1,param->powern);
  const double zeta2n = pow(zeta2,param->powern2);
  return param->zeta_shift_total*pow(pow(param->zeta_shift1+zeta1n, -1.0*param->powerp)+pow(param->zeta_shift2+zeta2n, -1.0*param->powerp), 1.0/(2.0*param->powern_del*param->powerp));
}


void PairTersoff_k2::ters_bij_d(const double zeta1, const double zeta2, Param *param, double& db1, double& db2)
{
  if (zeta1 > param->ca1 && zeta2 > param->ca1_2){
    if (zeta1 > zeta2){
      db1 = 0.0;
      db2 = -1.0*param->zeta_shift_total*param->powern2/(2.0*param->powern_del)*pow(zeta2, -1.0-param->powern2/(2.0*param->powern_del));
    }else{
      db1 = -1.0*param->zeta_shift_total*param->powern/(2.0*param->powern_del)*pow(zeta1, -1.0-param->powern/(2.0*param->powern_del));
      db2 = 0.0;
    }
  }
  else if (zeta1 < param->ca4 || zeta2 < param->ca4){
		db1 = 0.0;
    db2 = 0.0;
	}
  else{		  
    const double zeta1n = pow(zeta1,param->powern);
    const double zeta2n = pow(zeta2,param->powern2);
    const double zeta_d = pow(param->zeta_shift1+zeta1n, -1.0*param->powerp)+pow(param->zeta_shift2+zeta2n, -1.0*param->powerp);
    const double tmp = -0.5*param->zeta_shift_total/param->powern_del*pow(zeta_d, -1+1.0/(2.0*param->powern_del*param->powerp));
    db1 = tmp*param->powern*pow(param->zeta_shift1+zeta1n, -1.0-param->powerp)*zeta1n/zeta1;
    db2 = tmp*param->powern2*pow(param->zeta_shift2+zeta2n, -1.0-param->powerp)*zeta2n/zeta2;
  }
}

/* ---------------------------------------------------------------------- */

void PairTersoff_k2::ters_zetaterm_d(const double prefactor,
                                  double *rij_hat, const double rij,
                                  double *rik_hat, const double rik,
                                  double *dri, double *drj, double *drk,
                                  Param* const param, const double powermint, const double alpha,
                                  const double h, const double c, const double d, const double bigre_diff)
{
  double gijk,gijk_d,ex_delr,ex_delr_d,fc,dfc,cos_theta,tmp;
  double dcosdri[3],dcosdrj[3],dcosdrk[3];

  fc = ters_fc(1, rik,param);
  dfc = ters_fc_d(1, rik,param);

	double rdiff = rij-rik-bigre_diff;
  if (powermint == 3) tmp = alpha*pow(rdiff,3.0);
  else tmp = alpha * rdiff;

  if (tmp > 69.0776) ex_delr = 1.e30;
  else if (tmp < -69.0776) ex_delr = 0.0;
  else ex_delr = exp(tmp);

  if (powermint == 3)
    ex_delr_d = 3.0*alpha * pow(rdiff,2.0)*ex_delr;
  else ex_delr_d = alpha * ex_delr;

  cos_theta = vec3_dot(rij_hat,rik_hat);
  gijk = ters_gijk(cos_theta,h, c, d);
  gijk_d = ters_gijk_d(cos_theta,h,d);
  costheta_d(rij_hat,rij,rik_hat,rik,dcosdri,dcosdrj,dcosdrk);

  // compute the derivative wrt Ri
  // dri = -dfc*gijk*ex_delr*rik_hat;
  // dri += fc*gijk_d*ex_delr*dcosdri;
  // dri += fc*gijk*ex_delr_d*(rik_hat - rij_hat);

  vec3_scale(-dfc*gijk*ex_delr,rik_hat,dri);
  vec3_scaleadd(fc*gijk_d*ex_delr,dcosdri,dri,dri);
  vec3_scaleadd(fc*gijk*ex_delr_d,rik_hat,dri,dri);
  vec3_scaleadd(-fc*gijk*ex_delr_d,rij_hat,dri,dri);
  vec3_scale(prefactor,dri,dri);

  // compute the derivative wrt Rj
  // drj = fc*gijk_d*ex_delr*dcosdrj;
  // drj += fc*gijk*ex_delr_d*rij_hat;

  vec3_scale(fc*gijk_d*ex_delr,dcosdrj,drj);
  vec3_scaleadd(fc*gijk*ex_delr_d,rij_hat,drj,drj);
  vec3_scale(prefactor,drj,drj);

  // compute the derivative wrt Rk
  // drk = dfc*gijk*ex_delr*rik_hat;
  // drk += fc*gijk_d*ex_delr*dcosdrk;
  // drk += -fc*gijk*ex_delr_d*rik_hat;

  vec3_scale(dfc*gijk*ex_delr,rik_hat,drk);
  vec3_scaleadd(fc*gijk_d*ex_delr,dcosdrk,drk,drk);
  vec3_scaleadd(-fc*gijk*ex_delr_d,rik_hat,drk,drk);
  vec3_scale(prefactor,drk,drk);
}

/* ---------------------------------------------------------------------- */

void PairTersoff_k2::costheta_d(double *rij_hat, double rij,
                             double *rik_hat, double rik,
                             double *dri, double *drj, double *drk)
{
  // first element is devative wrt Ri, second wrt Rj, third wrt Rk

  double cos_theta = vec3_dot(rij_hat,rik_hat);

  vec3_scaleadd(-cos_theta,rij_hat,rik_hat,drj);
  vec3_scale(1.0/rij,drj,drj);
  vec3_scaleadd(-cos_theta,rik_hat,rij_hat,drk);
  vec3_scale(1.0/rik,drk,drk);
  vec3_add(drj,drk,dri);
  vec3_scale(-1.0,dri,dri);
}

/* ---------------------------------------------------------------------- */

double PairTersoff_k2::get_chi_i(int i)
{
  return pre_qeq_i[i].chi;
}

double PairTersoff_k2::get_bigj_ii(int i)
{
  return pre_qeq_i[i].bigj_ii;
}

const sparse_matrix& PairTersoff_k2::get_H() const
{
  return H;
}

