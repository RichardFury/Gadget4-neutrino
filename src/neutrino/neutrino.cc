/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Rui Hu, while the main part of GADGET4 N-body/SPH code is 
 * \copyright	developed by Volker Springel. Copyright (C) 2014-2020 by 
 * \copyright	Volker Springel (vspringel@mpa-garching.mpg.de) and all
 * \copyright   contributing authors.
 *******************************************************************************/

/*! \neutrino.cc 
* 
*	\brief module for neutrino evolution
*/

#ifdef NEUTRINO

#include <gsl/gsl_rng.h>
#include <stdio.h>

#include "../neutrino/neutrino.h"
#include "../data/allvars.h"
#include "../data/mymalloc.h"
#include "../data/dtypes.h"
#include "../logs/logs.h"
#include "../logs/timer.h"
#include "../system/system.h"
#include "../time_integration/timestep.h"

/** \brief Compute the evolution of neutrino.
* 
*/

void nusfr::NuParamCalculate(void)
{ 
	All.H0 = All.Hubble * All.HubbleParam * 1000. / (1E6 * 3.0857E16);
	All.Rhocr = (cc * 3.0 * All.H0 * All.H0) / (8.0 * M_PI * Gr);
    All.unittrans = pow(kb, 4) / ((pow(hbar, 3)) * (pow(c, 3)) * 2 * M_PI * M_PI);
    Check_Neutrinos();
}

void nusfr::Check_Neutrinos(void) 
{
    switch(All.LeptonAsymmetry)
    {
        case 2:
           All.Xi_2 = 0;
           All.Xi_1 = 0;
           break;

        case 1: 
           All.Xi_2 = Cal_Xi2(All.Xi_3);
           All.Xi_1 = Cal_Xi1(All.Xi_2, All.Xi_3);
           break;        

        case 0: 
           All.Xi2 = All.Xi_3;
           All.Xi1 = All.Xi_3;
           break;

        default:
           break;
    }

    switch(All.MassHierarchy)
    {
        case 0: /* Normal */
          All.Mass_2 = sqrt(All.Mass_1 * All.Mass_1 + 7.59e-5);
          All.Mass_3 = sqrt(All.Mass_2 * All.Mass_2 + 2.32e-3);
          break;

        case 1: /* Inverted */
          All.Mass_2 = sqrt(All.Mass_1 * All.Mass_1 + 2.32e-3);
          All.Mass_3 = sqrt(All.Mass_2 * All.Mass_2 + 7.59e-5);
          break;

        case 2: /* Identical */
          All.Mass_2 = All.Mass_1;
          All.Mass_3 = All.Mass_1;
          break;

        case 3: /* for only one sterile/active neutrino */
          All.Mass_2 = 0;
          All.Mass_3 = All.Mass_1;
          All.Mass_1 = 0;
          break;

        default:
          break;
    }
    
    printf("Mass1: %f Mass2: %f Mass3: %f Xi_1: %f Xi_2: %f Xi_3: %f \n", All.Mass_1, All.Mass_2, All.Mass_3, All.Xi_1, All.Xi_2,
           All.Xi_3);

    All.NumDensity0 = numdens(0.);
    All.NumDensity3 = numdens(All.Xi_3);
#ifdef STERILE
    All.NumDensity4 = numdens(All.Xi_4);
#endif  // STERILE

    switch(All.ExpanOn)
    {
        case 2: /* for only one sterile/active neutrino */
#ifdef STERILE
          All.Omega_Nu0_Expansion = neutrino_integration(1.0, All.Mass_3, All.Xi_3) + neutrino_integration(1.0, All.Mass_4, All.Xi_4);
          All.Omega_Nu0_Frstr     = neutrino_integration(1.0, 0., 0.) * 2;
#else
          All.Omega_Nu0_Expansion = neutrino_integration(1.0, All.Mass_3, All.Xi_3);
          All.Omega_Nu0_Frstr     = neutrino_integration(1.0, 0., 0.);
#endif  // STERILE
          All.NumDensity1 = 0;
          All.NumDensity2 = 0;
          break;

        case 1:
#ifdef STERILE
          All.Omega_Nu0_Expansion = neutrino_integration(1.0, All.Mass_4, All.Xi_4) + neutrino_integration(1.0, All.Mass_3, All.Xi_3) +
                                    neutrino_integration(1.0, All.Mass_2, All.Xi_2) + neutrino_integration(1.0, All.Mass_1, All.Xi_1);
          All.Omega_Nu0_Frstr     = neutrino_integration(1.0, 0., 0.) * 4;
#else
          All.Omega_Nu0_Expansion = neutrino_integration(1.0, All.Mass_3, All.Xi_3) + neutrino_integration(1.0, All.Mass_2, All.Xi_2) +
                                    neutrino_integration(1.0, All.Mass_1, All.Xi_1);
          All.Omega_Nu0_Frstr     = neutrino_integration(1.0, 0., 0.) * 3;
#endif  // STERILE
          All.NumDensity1         = numdens(All.Xi_1);
          All.NumDensity2         = numdens(All.Xi_2);
          break;

        case 0:
#ifdef STERILE
          All.Omega_Nu0_Expansion = neutrino_integration(1.0, 0., 0.) * 4;
          All.Omega_Nu0_Frstr     = neutrino_integration(1.0, 0., 0.) * 4;
#else
          All.Omega_Nu0_Expansion = neutrino_integration(1.0, 0., 0.) * 3;
          All.Omega_Nu0_Frstr     = neutrino_integration(1.0, 0., 0.) * 3;
#endif  // STERILE
          All.NumDensity1 = numdens(All.Xi_1);
          All.NumDensity2 = numdens(All.Xi_2);
          break;

        default:
          break;
    }
    
    switch(All.DeductionFromDE)
    {
        case 0:
          All.Omega2      = All.Omega0 - All.Omega_nu0_expan;
          All.OmegaLambda = All.OmegaLambda;
          break;

        case 1:
          All.Omega2      = All.Omega0;
          All.OmegaLambda -= All.Omega_nu0_expan;
          All.Omega0      += All.Omega_nu0_expan;
    }

    printf("Omega0 = %f Omega2 = %f Omega_Nu_Expansion = %f Omega_Nu_Frstr %f\t", All.Omega0, All.Omega2, All.Omega_Nu0_Expansion,
           All.Omega_Nu0_Frstr);
    printf("Omega_lambda = %f\n", All.OmegaLambda);

    printf("NumDensity %f %f %f %f\n", All.NumDensity0, All.NumDensity1, All.NumDensity2, All.NumDensity3);
}

double nusfr::hubble_function_nu(double a) 
{
  double h;
  double rhoneu;
  switch(All.expan_on)
    {
      case 2:
#ifdef STERILE
        rhoneu = neutrino_integration(a, All.Mass_3, All.Xi_3) + neutrino_integration(a, All.Mass_4, All.Xi_4);
#else
        rhoneu = neutrino_integration(a, All.Mass_3, All.Xi_3);
#endif  // STERILE
        break;

      case 1:
#ifdef STERILE
        rhoneu = neutrino_integration(a, All.Mass_1, All.Xi_1) + neutrino_integration(a, All.Mass_2, All.Xi_2) +
                 neutrino_integration(a, All.Mass_3, All.Xi_3) + neutrino_integration(a, All.Mass_4, All.Xi_4);
#else
        rhoneu = neutrino_integration(a, All.Mass_1, All.Xi_1) + neutrino_integration(a, All.Mass_2, All.Xi_2) +
                 neutrino_integration(a, All.Mass_3, All.Xi_3);
#endif  // STERILE
        break;

      case 0:
#ifdef STERILE
        rhoneu = neutrino_integration(a, 0., 0.) * 4;
#else
        rhoneu = neutrino_integration(a, 0., 0.) * 3;
#endif  // STERILE
        break;

      default:
        break;
    }

  h = All.Omega2 / (a * a * a) + (1 - All.Omega2 - All.OmegaLambda - All.Omega_Nu0_Expansion) / (a * a) + All.OmegaLambda + rhoneu;
  h = All.Hubble * sqrt(h);
  return h;
}