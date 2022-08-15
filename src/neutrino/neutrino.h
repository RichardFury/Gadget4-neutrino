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

#ifndef NEUTRINO_H
#define NEUTRINO_H



#include <fftw3.h>

typedef ptrdiff_t fft_ptrdiff_t;

#ifdef DOUBLEPRECISION_FFTW
typedef double fft_real;
typedef fftw_complex fft_complex;
#else
typedef float fft_real;
typedef fftwf_complex fft_complex;
#endif

#include "../data/simparticles.h"
#include "../mpi_utils/setcomm.h"
#include "../pm/pm_mpi_fft.h"

/*	Some constants used in the calculation of neutrino
*	However, this part is needed to be modified for consistence
*/
#define	hbar		      1.05457173E-34  // m^2*kg*s^-1
#define kb			      1.3806488E-23     // m^2*kg*s^-2*k^-1
#define Gr			      6.67384E-11       // m^3*kg^-1*s^-2
#define clight        2.9979e8           // m*s^-1
#define cclight	      8.9874e16
#define ktoev		      8.617e-5
#define mpc_to_m	    3.0857e22
#define ev_to_kg	    1.78266191E-36
#define mass_of_sun   1.989E30
#define relativistic  1

class nusfr : public pm_mpi_fft
{
 public:
  nusfr(MPI_Comm comm) : setcomm(comm), pm_mpi_fft(comm) {}
  
  void NeutrinoPkInit(void);

  void NeutrinoPkGenerate(fft_plan &my_plan, fft_complex *&fft_of_rhogrid);

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


  double phi(double q, double mass, double xi, double a);
  double phi_expansion(double q, double mass, double xi, double a);
  double phi_fit(double q, double mass, double a);
  
  double hubble(double a);
};

struct neu_params
{
  double x, y, z;
};

struct my_f_params
{
  double x, y, z;
  double w;
};

double Cal_Xi2(double xi3);

double Cal_Xi1(double xi2, double xi3);

double hubble_function_nu(double a);

double neutrino_integration(double a, double m, double xi);

double neutrino_partition(double pot, void* par);

void NuParamCalculate(void);

void Check_Neutrinos(void);

double numdens(double xi);

double numdens_integrand(double xi, void *par);

double phi_integrand(double u, void *par);

double phi_to_deduct(double u, void *par);

#endif  // NEUTRINO_H
