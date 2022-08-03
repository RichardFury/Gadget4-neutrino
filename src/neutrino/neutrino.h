/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Rui Hu, while the main part of GADGET4 N-body/SPH code is
 * \copyright	developed by Volker Springel. Copyright (C) 2014-2020 by
 * \copyright	Volker Springel (vspringel@mpa-garching.mpg.de) and all
 * \copyright   contributing authors.
 *******************************************************************************/

/*! \neutrino.cc
 *	\frstr.cc
 *
 *	\brief defines a class for dealing with neutrino evolution
 *  \brief defines some functions for calculating free streaming
 */

#ifdef NEUTRINO_H
#define NEUTRINO_H

#ifdef NEUTRINO

#include "../data/simparticles.h"
#include "../mpi_utils/setcomm.h"

/*	Some constants used in the calculation of neutrino
*	However, this part is needed to be modified for consistence
*/
#define	hbar		1.05457173E-34  // m^2*kg*s^-1
#define kb			1.3806488E-23     // m^2*kg*s^-2*k^-1
#define Gr			6.67384E-11       // m^3*kg^-1*s^-2
#define c			2.9979e8           // m*s^-1
#define cc			8.9874e16
#define ktoev		8.617e-5
#define mpc_to_m	3.0857e22
#define ev_to_kg	1.78266191E-36
#define mass_of_sun 1.989E30
#define relativistic 1

class nusfr : public setcomm
{
 public:
  nustr(MPI_Comm comm) : stecomm(comm) {}

  void NuParamCalculate(void);
  
  double cal_xi2(double xi3);

  double cal_xi1(double xi2, double xi3);

  double frstr(double k, double pk0_nu, double pk0_cdm, double pk1_cdm, double *array1, double *array2, double frstr_a0,
               double frstr_a1, double mass, double xi);

  double numdens(double xi);

  double a_to_s(int i, double *array1, double frstr_a0, double frstr_a1);

 private:

#define GRIDX ((PMGRID / LONG_X) * DBX + DBX_EXTRA)
#define GRIDY ((PMGRID / LONG_Y) * DBY + DBY_EXTRA)
#define GRIDZ ((PMGRID / LONG_Z) * DBZ + DBZ_EXTRA)

#if(GRIDX > 1024) || (GRIDY > 1024) || (GRIDZ > 1024)
  typedef long long large_array_offset; /* use a larger data type in this case so that we can always address all cells of the 3D grid
                                           with a single index */
#else
  typedef int large_array_offset;
#endif

  double hubble_function_nu(double a);

  struct my_f_params
  {
    double x, y, z;
    double w;
  };

  struct neu_params
  {
    double x, y, z;
  };

  double phi_integrand(double u, void *par);

  double phi(double q, double mass, double xi, double a);
  double phi_expansion(double q, double mass, double xi, double a);
  double phi_fit(double q, double mass, double a);
  double phi_to_deduct(double u, void *par);
  double hubble(double a);
  double neutrino_integration(double a, double m, double xi);

  
  double numdens_integrand(double xi, void *par);


};

#endif  // NEUTRINO

#endif  // NEUTRINO_H
