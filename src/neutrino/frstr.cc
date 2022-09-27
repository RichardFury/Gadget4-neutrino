/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Rui Hu, while the main part of GADGET4 N-body/SPH code is
 * \copyright	developed by Volker Springel. Copyright (C) 2014-2020 by
 * \copyright	Volker Springel (vspringel@mpa-garching.mpg.de) and all
 * \copyright   contributing authors.
 *******************************************************************************/

/*! \frstr.cc
 *
 *	\brief function for free streaming calculations and xi calculation.
 */


#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <string.h>
#include <gsl/gsl_integration.h>

#include "../data/allvars.h"
#include "../neutrino/neutrino.h"

double Cal_Xi2(double xi3)
{
  double AA, BB, CC, DD, t12, t13, s13, s23, c23, s12, c12, c13, r23;
  double xi2, L3, L2;
  int i;

  s12 = sqrt(0.304);
  s23 = sqrt(0.51);
  s13 = sqrt(0.0219);
  c23 = sqrt(1. - 0.51);
  c12 = sqrt(1. - 0.304);
  c13 = sqrt(1. - 0.0219);
  t12 = s12 / c12;
  t13 = s13 / c13;

  AA = c23 * ((1. - t12 * t12) * c23 - 2. * s13 * s23 * t12);
  BB = ((1. - t13 * t13) * s23 * s23 - t12 * t13 * t13 * c23 * (2. * s13 * s23 + t12 * c23));
  CC = s23 * ((1. - t12 * t12) * s23 + 2. * s13 * c23 * t12);
  DD = ((1. - t13 * t13) * c23 * c23 + t12 * t13 * t13 * s23 * (2. * s13 * c23 - t12 * s23));

  r23 = (DD - BB) / (AA - CC);

  L3 = xi3 * (xi3 * xi3 + M_PI * M_PI);
  L2 = r23 * L3;

  xi2 = L2 / (M_PI * M_PI);
  for(i = 0; i < 10; i++)
    {
      xi2 = L2 / (xi2 * xi2 + M_PI * M_PI);
      //   printf("xi2 %f\n", xi2);
    }

  return xi2;
}

double Cal_Xi1(double xi2, double xi3) {
  double s13, s12, s23, c12, c13, c23;
  double xi1, L3, L2, L1;
  int i;

  s12 = sqrt(0.304);
  s23 = sqrt(0.51);
  s13 = sqrt(0.0219);
  c23 = sqrt(1. - 0.51);
  c12 = sqrt(1. - 0.304);
  c13 = sqrt(1. - 0.0219);

  L3 = xi3 * (xi3 * xi3 + M_PI * M_PI);
  L2 = xi2 * (xi2 * xi2 + M_PI * M_PI);

  L1  = -(s13 * s13 * L3 + c13 * c13 * s12 * s12 * L2) / (c13 * c13 * c12 * c12);
  xi1 = L1 / (M_PI * M_PI);

  for(i = 0; i < 10; i++)
    {
      xi1 = L1 / (xi1 * xi1 + M_PI * M_PI);
      //  printf("xi1 %f\n", xi1);
    }

  return xi1;
}

double neutrino_partition(double pot, void* par) {
  struct neu_params* params = (struct neu_params*)par;
  double part;
  double Tneu    = (params->x);
  double mass_nu = (params->y);
  double nu_xi   = (params->z);
  double mass_n;

  mass_n = mass_nu / ktoev;
  part   = pow(pot, 3) * sqrt(1 + pow(mass_n / (pot * Tneu), 2)) / (exp((pot)-nu_xi) + 1.0) +
         pow(pot, 3) * sqrt(1 + pow(mass_n / (pot * Tneu), 2)) / (exp((pot) + nu_xi) + 1.0);
  // printf("massn %f mass nu %f, ktoev %f part %f Tneu %f\n", mass_n, mass_nu, ktoev, part,Tneu);

  return part;
}

double neutrino_integration(double a, double m, double xi)
{
#define WORKSIZE2 100000
  gsl_function F2;
  gsl_integration_workspace *workspace2;

  double inte_result, inte_abserr;
  double roneu;
  double Tneu;
  Tneu                    = All.T_Neutrino_0 / a;
  struct neu_params alpha = {Tneu, m, xi};

  workspace2  = gsl_integration_workspace_alloc(WORKSIZE2);
  F2.function = &neutrino_partition;
  F2.params   = &alpha;

  gsl_integration_qagiu(&F2, 0.0, 0, 1.0e-8, WORKSIZE2, workspace2, &inte_result, &inte_abserr);
  roneu = inte_result * All.Unittrans * pow(Tneu, 4);
  roneu = roneu / All.Rhocr;

   //printf("mass %f xi %f a %f roneu %f rocr %f unittrans %f inte_result %f\n", m, xi, a, roneu, All.Rhocr *1e15, All.Unittrans*1e15,
   //inte_result);

  gsl_integration_workspace_free(workspace2);
  return roneu;
}

double nusfr::hubble(double a)
{
  double hub, rrp, rom, ror, rolambda, roneu, rototal;
  int i;
  rom      = All.Omega2 * All.Rhocr / (a * a * a);
  rolambda = All.OmegaLambda * All.Rhocr;
  if(All.ExpanOn == 1)
    {
      roneu = 0.0;
      for(i = 0; i < All.NNeutrino; i++)
        {
#ifdef STERILE
          if(i == STERILE)
            {
              roneu += All.Neff * neutrino_integration(a, All.NuMass[i], All.Xi[i]);
            }
#else
          roneu += neutrino_integration(a, All.NuMass[i], All.Xi[i]);
#endif
        }
    }

  if(All.ExpanOn == 0)
    {
      roneu = neutrino_integration(a, 0., 0.) * All.NNeutrino;
    }

  roneu *= All.Rhocr;

  rototal = rom + roneu + rolambda;
  hub     = sqrt(8 * M_PI * Gr * rototal / 3) / clight;

  return hub;
}

double nusfr::a_to_s(int i, double *array1, double frstr_a0, double frstr_a1)
{
  int j;
  double temp;
  double a_spacing = (frstr_a1 - frstr_a0) / (double)(All.FrstrInterval);
  temp             = 0.;
  /*
    if(i == 1){
    temp = 0.5 / pow(array1[i-1], 3) / hubble(array1[i-1]) + 0.5 / pow(array1[i], 3) / hubble(array1[i]);
  }

  else{
    if(i%2 == 0){
      for(j=1;j<i;j++){
            if(j%2 == 0){
                temp = temp + 2. / pow(array1[j], 3) / hubble(array1[j]);
              }
            else{
                temp = temp + 4. / pow(array1[j], 3) / hubble(array1[j]);
            }
          }
      temp = temp + 1. / pow(array1[i], 3) / hubble(array1[i]);
      temp = temp / 3.;
    }

    else{
      for(j=1;j<i-1;j++){
            if(j%2 == 0){
                temp = temp + 2. / pow(array1[j], 3) / hubble(array1[j]);
              }
            else{
                temp = temp + 4. / pow(array1[j], 3) / hubble(array1[j]);
            }
          }
      temp = temp + 1. / pow(array1[i-1], 3) / hubble(array1[i-1]);
      temp = temp / 3.;
      temp = temp + 0.5 / pow(array1[i-1], 3) / hubble(array1[i-1]) + 0.5 / pow(array1[i], 3) / hubble(array1[i]);
    }
  }

   //printf("hubble %f temp %f aspacing %f a %f\n", hubble(array1[i-1])*1e18, temp, a_spacing, array1[i]);
  */

  for(j = 1; j <= i; j++)
    {
      temp += 0.5 / pow(array1[j - 1], 3) / hubble(array1[j - 1]) + 0.5 / pow(array1[j], 3) / hubble(array1[j]);
    }

  return temp * a_spacing;
}

double nusfr::frstr(double k, double pk0_nu, double pk0_cdm, double pk1_cdm, double *array1, double *array2, double frstr_a0, double frstr_a1,
             double mass, double xi)
{
  double *phi2;
  int i, knum, j, m;
  int Volterra_iteration_num = 1;
  double a_spacing           = (frstr_a1 - frstr_a0) / (double)(All.FrstrInterval);

  double phi0, phi1;
  phi2 = (double *)malloc((All.FrstrInterval + 1) * sizeof(double));

  double *k1;
  k1 = (double *)malloc((All.FrstrInterval + 1) * sizeof(double));

  double *k2;
  k2 = (double *)malloc((All.FrstrInterval + 1) * sizeof(double));

  int integration_equation_solver = 1;

  for(i = 0; i < All.FrstrInterval; i++)
    {
      // printf("arr2 %lf arr1 %lf k %f array 2 %f\n", k*(array2[All.frstr_interval] - array2[i]), array1[i], k, array2[i]);
      if(All.PhiParam == 1)
        {
          // phi1[i] = phi(k*(array2[All.frstr_interval] - array2[i]), array1[i]) * mpc_to_m / pow(c, 3);
          phi2[i] = phi(k * array2[i], mass, xi, array1[i]);
        }

      if(All.PhiParam == 2)
        {
          // phi1[i] = phi_expansion(k*(array2[All.frstr_interval] - array2[i]), array1[i]) * mpc_to_m / pow(c, 3);
          phi2[i] = phi_expansion(k * array2[i], mass, xi, array1[i]);
        }

      if(All.PhiParam == 3)
        {
          // phi1[i] = phi_fit(k*(array2[All.frstr_interval] - array2[i]), array1[i]) * mpc_to_m / pow(c, 3);
          phi2[i] = phi_fit(k * array2[i], mass, array1[i]);
        }

      // printf("phi2 %.10f i %d q %f phi1*Gr*All.rocr%.15f\n", phi2[i], i, array2[i], phi1[i] * Gr * All.rocr);

      k1[i] = 0.;
      k2[i] = 0.;
    }

  // printf("phi loop ended\n");
  if(All.PhiParam == 1)
    {
      phi2[All.FrstrInterval] = phi(k * array2[All.FrstrInterval], mass, xi, array1[All.FrstrInterval]);
    }

  if(All.PhiParam == 2)
    {
      phi2[All.FrstrInterval] = phi_expansion(k * array2[All.FrstrInterval], mass, xi, array1[All.FrstrInterval]);
    }

  if(All.PhiParam == 3)
    {
      phi2[All.FrstrInterval] = phi_fit(k * array2[All.FrstrInterval], mass, array1[All.FrstrInterval]);
    }

  // calculate 3d delta k
  if(integration_equation_solver == 1)
    {
      // iteration to solve the Volterra equation, for the discrete inte part, trapazoidal rule is temporarily used here
      double f0, f1, m0, m1, f_temp, m_temp, m_temp2, f_temp2;

      k1[0]        = 0.;
      double munit = 4 * M_PI * Gr * All.Rhocr * a_spacing;
      for(i = 1; i <= All.FrstrInterval; i++)
        {
          m_temp = 0.;
          for(j = 1; j <= i; j++)
            {
              double temp_omegam;

              temp_omegam = All.Omega0 - All.Omega_Nu0_Frstr;

              if(All.PhiParam == 1)
                {
                  phi0 = phi(k * (array2[i] - array2[j - 1]), mass, xi, array1[j - 1]) * mpc_to_m / pow(clight, 3);
                  phi1 = phi(k * (array2[i] - array2[j]), mass, xi, array1[j]) * mpc_to_m / pow(clight, 3);
                }

              if(All.PhiParam == 2)
                {
                  phi0 = phi_expansion(k * (array2[i] - array2[j - 1]), mass, xi, array1[j - 1]) * mpc_to_m / pow(clight, 3);
                  phi1 = phi_expansion(k * (array2[i] - array2[j]), mass, xi, array1[j]) * mpc_to_m / pow(clight, 3);
                  // printf("phi0 %f i %d j %d k %f\n", phi0, i, j, k);
                }

              if(All.PhiParam == 3)
                {
                  phi0 = phi_fit(k * (array2[i] - array2[j - 1]), mass, array1[j - 1]) * mpc_to_m / pow(clight, 3);
                  phi1 = phi_fit(k * (array2[i] - array2[j]), mass, array1[j]) * mpc_to_m / pow(clight, 3);
                  // printf("phi0 %f i %d j %d ti %f tj %f\n", phi0, i, j, array2[i], array2[j-1]);
                }

              m0 = phi0 * array1[j - 1] * (array2[i] - array2[j - 1]) * temp_omegam *
                   ((pk1_cdm - pk0_cdm) * (j - 1) / All.FrstrInterval + pk0_cdm) / (pow(array1[j - 1], 3) * hubble(array1[j - 1]));

              m1 = phi1 * array1[j] * (array2[i] - array2[j]) * temp_omegam *
                   ((pk1_cdm - pk0_cdm) * j / All.FrstrInterval + pk0_cdm) / (pow(array1[j], 3) * hubble(array1[j]));
              // printf("phi0 %f phi1 %f\n", phi0, phi1);
              m_temp = m_temp + 0.5 * (m0 + m1);
            }
          m_temp2 = m_temp * munit;

          k1[i] = m_temp2 + phi2[i] * pk0_nu;
          // printf("m_temp2 %.10f, phi2[i] %f q %f\n", m_temp2, phi2[i], k*array2[i]);
        }
      // k0 is the higher redshift delta in k space at (kx, ky, kz)

      for(i = 0; i <= All.FrstrInterval; i++)
        {
          k2[i] = k1[i];
        }

      for(m = 0; m < Volterra_iteration_num; m++)
        {
          double onu_temp0, onu_temp1;
          for(i = 1; i <= All.FrstrInterval; i++)
            {
              f_temp = 0.;
              for(j = 1; j <= i; j++)
                {
                  onu_temp0 = neutrino_integration(array1[j - 1], mass, xi);
                  onu_temp1 = neutrino_integration(array1[j], mass, xi);

                  onu_temp0 *= pow(
                      array1[j - 1],
                      3);  // this is because the neutrino_integration func returns Omega in unit of rocr at 0 in terms of percentage
                  onu_temp1 *= pow(array1[j], 3);

                  if(All.PhiParam == 1)
                    {
                      phi0 = phi(k * (array2[i] - array2[j - 1]), mass, xi, array1[j - 1]) * mpc_to_m / pow(clight, 3);
                      phi1 = phi(k * (array2[i] - array2[j]), mass, xi, array1[j]) * mpc_to_m / pow(clight, 3);
                    }

                  if(All.PhiParam == 2)
                    {
                      phi0 = phi_expansion(k * (array2[i] - array2[j - 1]), mass, xi, array1[j - 1]) * mpc_to_m / pow(clight, 3);
                      phi1 = phi_expansion(k * (array2[i] - array2[j]), mass, xi, array1[j]) * mpc_to_m / pow(clight, 3);
                    }

                  if(All.PhiParam == 3)
                    {
                      phi0 = phi_fit(k * (array2[i] - array2[j - 1]), mass, array1[j - 1]) * mpc_to_m / pow(clight, 3);
                      phi1 = phi_fit(k * (array2[i] - array2[j]), mass, array1[j]) * mpc_to_m / pow(clight, 3);
                    }

                  f0 = phi0 * array1[j - 1] * (array2[i] - array2[j - 1]) * onu_temp0 * k2[j - 1] /
                       (pow(array1[j - 1], 3) * hubble(array1[j - 1]));
                  f1 = phi1 * array1[j] * (array2[i] - array2[j]) * onu_temp1 * k2[j] / (pow(array1[j], 3) * hubble(array1[j]));

                  f_temp = f_temp + 0.5 * (f0 + f1);
                }
              f_temp2 = f_temp * munit;

              k2[i] = f_temp2 + k1[i];
              // printf("ftemp %lf k2%lf k1%lf f0%lf\n", f_temp2, k2[i], k1[i], f0);
            }
          // printf("Volterra eqn no %d, delta k nv %.8f k1 %.8f knv0 %lf k1part %.8f\n", m, k2[All.frstr_interval],
          // k1[All.frstr_interval], pk0_nu, phi2[All.frstr_interval] * pk0_nu);
        }
    }
  return k2[All.FrstrInterval];
}

double phi_integrand(double u, void *par)
{
  double up, down1, down2, T;

  int unitless_inte          = 1;
  struct my_f_params *params = (struct my_f_params *)par;

  double a    = (params->x);
  double q    = (params->y);
  double mass = (params->z);
  double xi   = (params->w);

  T             = All.T_Neutrino_0 / a;
  double con_ul = T * ktoev / mass;

  if(unitless_inte == 1)
    {
      if(relativistic == 0)
        {
          up    = u * (sin(u * q * con_ul) + 1.);
          down1 = exp(u - xi) + 1.;
          down2 = exp(u + xi) + 1.;
        }
      if(relativistic == 1)
        {
          up    = u / sqrt(1. - u * u * con_ul * con_ul) * (sin(u / sqrt(1. - u * u * con_ul * con_ul) * q * con_ul) + 1.);
          down1 = exp(u / sqrt(1. - u * u * con_ul * con_ul) - xi) + 1.;
          down2 = exp(u / sqrt(1. - u * u * con_ul * con_ul) + xi) + 1.;
        }
      // up = u * sin(u * q * con_ul);
      // expand in taylor
    }

  return up / down1 ;
}

double phi_to_deduct(double u, void *par)
{
  double up, down1, down2, T;

  struct my_f_params *params = (struct my_f_params *)par;

  double a    = (params->x);
  double q    = (params->y);
  double mass = (params->z);
  double xi   = (params->w);

  T             = All.T_Neutrino_0 / a;
  double con_ul = T * ktoev / mass;
  // up = u * sin(u * q * con_ul);
  if(relativistic == 0)
    {
      up    = u;
      down1 = exp(u - xi) + 1.;
      down2 = exp(u + xi) + 1.;
    }
  if(relativistic == 1)
    {
      up    = u / sqrt(1. - u * u * con_ul * con_ul);
      down1 = exp(u / sqrt(1. - u * u * con_ul * con_ul) - xi) + 1.;
      down2 = exp(u / sqrt(1. - u * u * con_ul * con_ul) + xi) + 1.;
    }

  return up / down1 ;
}


double nusfr::phi(double q, double mass, double xi, double a)
{
  double inte_result, inte_abserr;
  double inte_result2, inte_abserr2;
  double normal_factor;
  double TT = All.T_Neutrino_0 / a;
#define WORKSIZE2 100000
  gsl_function F2;
  gsl_function F3;
  gsl_function F4;
  gsl_integration_workspace *workspace2;
  gsl_integration_workspace *workspace3;
  gsl_integration_workspace *workspace4;

  struct my_f_params alpha = {a, q, mass, xi};

  workspace2  = gsl_integration_workspace_alloc(WORKSIZE2);
  workspace3  = gsl_integration_workspace_alloc(WORKSIZE2);
  workspace4  = gsl_integration_workspace_alloc(WORKSIZE2);
  F2.function = &phi_integrand;
  F2.params   = &alpha;

  F3.function = &phi_to_deduct;
  F3.params   = &alpha;

  // F4.function = &normal_integrand;
  // F4.params = &TT;

  double con_ul = (All.T_Neutrino_0 / a) * ktoev / mass;
  // gsl_integration_qagiu(&F2, 0.0, 0, 1.0e-2, WORKSIZE2, workspace2, &inte_result, &inte_abserr);
  gsl_integration_qag(&F2, 0.0, 2., 0, 1.0e-2, WORKSIZE2, GSL_INTEG_GAUSS61, workspace2, &inte_result, &inte_abserr);
  gsl_integration_qag(&F3, 0.0, 2., 0, 1.0e-2, WORKSIZE2, GSL_INTEG_GAUSS61, workspace3, &inte_result2, &inte_abserr2);
  // gsl_integration_qag(&F4, 0.0, 2., 0, 1.0e-2, WORKSIZE2, GSL_INTEG_GAUSS61, workspace4, &normal_factor, &inte_abserr2);
  // gsl_integration_qagiu(&F3, 0.0, 0, 1.0e-2, WORKSIZE2, workspace3, &inte_result2, &inte_abserr2);
  gsl_integration_workspace_free(workspace2);
  gsl_integration_workspace_free(workspace3);
  gsl_integration_workspace_free(workspace4);
  // printf("phi check point 2 %.20f %lf\n", inte_result, q);

  int unitless_inte = 1;
  inte_result       = inte_result - inte_result2;
  double zeta3      = 1.202056;

  if(unitless_inte == 1)
    {
      // inte_result = inte_result * 4. / 3. / q / con_ul / zeta3;
      inte_result = inte_result / q / con_ul / normal_factor;
    }
  // question why integrate to 2.?
  return inte_result;
}

double nusfr::phi_expansion(double q, double mass, double xi, double a)
{
  double Tneu  = All.T_Neutrino_0 / a;
  double zeta3 = 1.202056;
  int i, n;
  double sum, temp, x, sum2;
  double con_ul = Tneu * ktoev / mass;

  sum  = 0.;
  sum2 = 0.;
  x    = q * con_ul;
  n    = 500;

  if(q == 0.)
    {
      sum2 == 0.;
    }

  else
    {
      if(xi == 0.)
        {
          for(i = 1; i <= n; i++)
            {
              temp = pow(-1., i + 1) * i / pow((i * i + x * x), 2);
              sum += temp;
            }
          sum2 = sum;
          sum2 *= 4. / (3. * zeta3);
        }

      if(xi != 0.)
        {
          double temp_dens;

          for(i = 0; i < All.NNeutrino; i++)
            {
              if(xi == All.Xi[i])
                temp_dens = All.NumDensity[i];
            }
          xi = abs(xi);
          double c1, a1, a2, a3, a4, sum, a5, temp1, temp2, err;
          a1 = 0.;
          a2 = 0.;
          a3 = 0.;
          a4 = 0.;
          a5 = 0.;

          for(i = 1; i <= n; i++)
            {
              a1 += pow(-1., i + 1) * i / (i * i + x * x);
              a2 += pow(-1., i + 1) * x / (i * i + x * x);
              a3 += pow(-1., i + 1) * (i * i - x * x) / pow((i * i + x * x), 2);
              a4 += pow(-1., i + 1) * (2. * x * i) / pow((i * i + x * x), 2);
              a5 += pow(-1., i + 1) * 2. * i * x * exp(-i * xi) / pow((i * i + x * x), 2);
            }

          struct my_f_params alpha = {a, q, mass, xi};
#define WORKSIZE2 100000
          gsl_function F2;
          gsl_function F3;

          gsl_integration_workspace *workspace2;
          gsl_integration_workspace *workspace3;
          workspace2 = gsl_integration_workspace_alloc(WORKSIZE2);
          workspace3 = gsl_integration_workspace_alloc(WORKSIZE2);

          F2.function = &phi_integrand;
          F2.params   = &alpha;

          F3.function = &phi_to_deduct;
          F3.params   = &alpha;

          gsl_integration_qag(&F2, 0.0, xi, 0, 1.0e-6, WORKSIZE2, GSL_INTEG_GAUSS61, workspace2, &temp1, &err);
          gsl_integration_qag(&F3, 0.0, xi, 0, 1.0e-6, WORKSIZE2, GSL_INTEG_GAUSS61, workspace3, &temp2, &err);
          gsl_integration_workspace_free(workspace2);
          gsl_integration_workspace_free(workspace3);

          c1 = temp1 - temp2;

          /*int num_inte_step = 30;
          c1 = 0.;
          //double *num_inte_temp;
          double temp_x0, temp_x1, temp_func;
          //num_inte_temp = (double*) malloc((num_inte_step + 1) * sizeof(double));

          for(i=0;i<num_inte_step;i++){
              temp_x0 = xi * i / num_inte_step;
              temp_x1 = xi * (i+1) / num_inte_step;
              temp_func = 0.25 * (temp_x0 * sin(x * temp_x0) / (exp(temp_x0 - xi) + 1) + temp_x1 * sin(x * temp_x1) / (exp(temp_x1 -
          xi) + 1)); temp_func += 0.25 * (temp_x0 * sin(x * temp_x0) / (exp(temp_x0 + xi) + 1) + temp_x1 * sin(x * temp_x1) /
          (exp(temp_x1 + xi) + 1)); c1 += temp_func * xi / num_inte_step;
          }*/

          sum = c1 + xi * cos(x * xi) * a2 + xi * sin(x * xi) * a1 + cos(x * xi) * a4 + sin(x * xi) * a3 + a5;
          // printf("c1 %f a2 %f a1 %f a3 %f a4 %f a5 %f\n", c1 * 1e10, xi * cos(x*xi)*a2*1e10, xi * sin(x*xi)*a1*1e10,
          // sin(x*xi)*a3*1e10, cos(x*xi)*a4*1e10, a5*1e10);
          sum2 = sum / 2. / x / 2.;  // the second 2 comes from the nu-antinu pair, that in the normalization part cancelled
          // printf("sum %f x %f\n", sum, x);
          sum2 *= 4. / temp_dens;
        }
    }

  // sum2 *= 4. / (3. * zeta3) * All.numdens0 / All.numdens1;

  if(sum2 > 1.)
    {
      sum2 = 0.99999;
    }

  return sum2;
}

double nusfr::phi_fit(double q, double mass, double a)
{
  double Tneu   = All.T_Neutrino_0 / a;
  double zeta3  = 1.202056;
  double con_ul = Tneu * ktoev / mass;
  double fit, x;

  x = q * con_ul;

  fit = (1. + 0.0168 * x * x + 0.0407 * pow(x, 4)) / (1. + 2.1734 * x * x + 1.6787 * pow(x, 4.1811) + 0.1467 * pow(x, 8));

  return fit;
}

double numdens_integrand(double x, void *par)
{
  double xi = *(double *)par;
  return x * x / (exp(x - xi) + 1.) + x * x / (exp(x + xi) + 1.);
}

double numdens(double xi)
{
#define WORKSIZE2 100000
  gsl_function F;
  gsl_integration_workspace *workspace2;

  double inte_result, inte_abserr;

  workspace2 = gsl_integration_workspace_alloc(WORKSIZE2);
  F.function = &numdens_integrand;
  F.params   = &xi;

  gsl_integration_qagiu(&F, 0.0, 0, 1.0e-8, WORKSIZE2, workspace2, &inte_result, &inte_abserr);
  gsl_integration_workspace_free(workspace2);

  return inte_result;
}


