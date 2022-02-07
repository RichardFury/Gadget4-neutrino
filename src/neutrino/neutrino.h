/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Rui Hu, while the main part of GADGET4 N-body/SPH code is
 * \copyright	developed by Volker Springel. Copyright (C) 2014-2020 by
 * \copyright	Volker Springel (vspringel@mpa-garching.mpg.de) and all
 * \copyright   contributing authors.
 *******************************************************************************/

/*! \neutrino.cc
 *
 *	\brief defines a class for dealing with neutrino evolution
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

class nusfr : public setcomm
{
 public:
  nustr(MPI_Comm comm) : stecomm(comm) {}

  void NuParamCalculate(void);
  double hubble_function_nu(double a);

};

#endif  // NEUTRINO

#endif  // NEUTRINO_H
