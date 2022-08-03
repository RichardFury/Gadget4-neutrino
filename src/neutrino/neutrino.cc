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

#include "../neutrino/neutrino.h"

#include <gsl/gsl_rng.h>
#include <stdio.h>

#include "../data/allvars.h"
#include "../data/dtypes.h"
#include "../data/mymalloc.h"
#include "../logs/logs.h"
#include "../logs/timer.h"
#include "../system/system.h"
#include "../time_integration/timestep.h"

/** \brief Compute the evolution of neutrino.
 *
 */

/* This function aims to calculate the modified cosmological parameters, Since there are neutrinos that modified the standard
 * cosmological model.
 *
 */

#define GRIDz (GRIDZ / 2 + 1)
#define GRID2 (2 * GRIDz)

void nusfr::NuParamCalculate(void)
{
  All.H0        = All.Hubble * All.HubbleParam * 1000. / (1E6 * 3.0857E16);
  All.Rhocr     = (cc * 3.0 * All.H0 * All.H0) / (8.0 * M_PI * Gr);
  All.Unittrans = pow(kb, 4) / ((pow(hbar, 3)) * (pow(c, 3)) * 2 * M_PI * M_PI);
  Check_Neutrinos();
}

void nusfr::Check_Neutrinos(void)
{
  switch(All.LeptonAsymmetry)
    {
      case 2:
        All.Xi[1] = 0;
        All.Xi[0] = 0;
        break;

      case 1:
        All.Xi[1] = Cal_Xi2(All.Xi[2]);
        All.Xi[0] = Cal_Xi1(All.Xi[1], All.Xi[2]);
        break;

      case 0:
        All.Xi[1] = All.Xi[2];
        All.Xi[0] = All.Xi[2];
        break;

      default:
        break;
    }

  switch(All.MassHierarchy)
    {
      case 0: /* Normal */
        All.NuMass[1] = sqrt(All.NuMass[0] * All.NuMass[0] + 7.59e-5);
        All.NuMass[2] = sqrt(All.NuMass[1] * All.NuMass[1] + 2.32e-3);
        break;

      case 1: /* Inverted */
        All.Mass_2 = sqrt(All.Mass_1 * All.Mass_1 + 2.32e-3);
        All.Mass_3 = sqrt(All.Mass_2 * All.Mass_2 + 7.59e-5);
        break;

      case 2: /* Identical */
        All.NuMass[1] = All.NuMass[0];
        All.NuMass[2] = All.NuMass[0];
        break;

      default:
        break;
    }

  for(int i = 0; i < All.NNeutrino; i++)
    {
      printf("Mass%d: %f\t Xi%d: %f\n", i + 1, All.NuMass[i], i + 1, All.Xi[i]);
      All.NumDensity[i] = numdens(All.Xi[i]);
    }

  All.NumDensity_0 = numdens(0.);

  if(All.ExpanOn)
    {
      for(int i = 0; i < All.NNeutrino; i++)
        {
          All.Omega_Nu0_Expansion += neutrino_integration(1.0, All.NuMass[i], All.Xi[i]);
        }
    }
  else
    {
      All.Omega_Nu0_Expansion = neutrino_integration(1.0, 0., 0.) * All.NNeutrino;
    }
  All.Omega_Nu0_Frstr = neutrino_integration(1.0, 0., 0.) * All.NNeutrino;

  if(All.DeductionFromDE)
    {
      All.Omega2 = All.Omega0;
      All.OmegaLambda -= All.Omega_nu0_expan;
      All.Omega0 += All.Omega_nu0_expan;
    }
  else
    {
      All.Omega2      = All.Omega0 - All.Omega_nu0_expan;
      All.OmegaLambda = All.OmegaLambda;
    }

  printf("Omega0 = %f Omega2 = %f Omega_Nu_Expansion = %f Omega_Nu_Frstr %f\t", All.Omega0, All.Omega2, All.Omega_Nu0_Expansion,
         All.Omega_Nu0_Frstr);
  printf("Omega_lambda = %f\n", All.OmegaLambda);

  printf("NumDensity %f", All.NumDensity_0);
  for(int i = 0; i < All.NNeutrino; i++)
    {
      printf("\t%f", All.NumDensity[i]);
    }
  printf("\n");
}


/* This function aims to calculate the hubble parameter at the current timestep.
 */
double nusfr::hubble_function_nu(double a)
{
  double h;
  double rhoneu;
  if(All.expan_on)
    {
      rhoneu = 0.;
      for(int i = 0; i < All.NNeutrino; i++)
        {
          rhoneu += neutrino_integration(a, All.NuMass[i], All.Xi[i]);
        }
    }
  else
    {
      rhoneu = neutrino_integration(a, 0., 0.) * All.NNeutrino;
    }
  h = All.Omega2 / (a * a * a) + (1 - All.Omega2 - All.OmegaLambda - All.Omega_Nu0_Expansion) / (a * a) + All.OmegaLambda + rhoneu;
  h = All.Hubble * sqrt(h);
  return h;
}

void nusfr::NeutrinoPkInit(void)
{
  double k_interval = 1.005;
  int num_kbins = (int)(log(sqrt(3.) * PMGRID / 0.95) / log(k_interval));
  if(All.NumCurrentTiStep == 0)
    {
      old_pk_b    = (double *)malloc(num_kbins * sizeof(double));
      old_pk_nu_b = (double *)malloc(num_kbins * All.NNeutrino * sizeof(double));
      count_b     = (double *)malloc(num_kbins * sizeof(double));
      k_array0    = (double *)malloc((num_kbins + 1) * sizeof(double));
      k_array     = (double *)malloc(num_kbins * sizeof(double));

      int ii;
      double ratio_temp;
      int rd, rd_size;

      File *fd;
      double kkk, ppp;

      if(!(fd = fopen(All.Ratio_Nu_CDM_Txt, "r")))
        {
          printf("can't read input spectrum in file '%s'\n", All.ratio_nu_cdm_txt);
        }

      rd_size = 0;
      do
        {
          if(fscanf(fd, " %lg %lg", &kkk, &ppp) == 2)
            {
              rd_size++;
            }
          else
            {
              break;
            }
        }
      while(1);
      fclose(fd);

      rd_array_k  = (double *)malloc(rd_size * sizeof(double));
      rd_array_pk = (double *)malloc(rd_size * sizeof(double));

      rd_size = 0;
      do
        {
          if(fscanf(fd, " %lg %lg", &kkk, &ppp) == 2)
            {
              rd_array_k[rd_size]  = kkk;
              rd_array_pk[rd_size] = ppp;
              rd_size++;
            }
          else
            {
              break;
            }
        }
      while(1);
      fclose(fd);

      FILE *outtimetxt;
      double aaa;
      int output_no;

      if(!(fd = fopen(All.OutputListFilename, "r")))
        {
          printf("can't read input spectrum in file '%s'\n", All.OutputListFilename);
        }

      do
        {
          if(fscanf(fd, "%lg", &aaa) == 1)
            {
              output_no++;
            }
          else
            break;
        }
      while(1);

      output_time_array = (double *)malloc(output_no * sizeof(double));

      output_no = 0;

      if(!(fd = fopen(All.OutputListFilename, "r")))
        {
          printf("can't read input spectrum in file '%s'\n", All.OutputListFilename);
        }

      do
        {
          if(fscanf(fd, "%lg", &aaa) == 1)
            {
              output_time_array[output_no] = aaa;
              output_no++;
            }
          else
            break;
        }
      for(ii = 0; ii < num_kbins; ii++)
        {
          old_pk_b[ii] = 0.;
          for(i = 0; i < All.NNeutrino; i++)
            {
              old_pk_nu_b[ii * All.NNeutrino + i] = 0.;
            }
        }
    }
}

void nusfr::NeutrinoPkGenerate(fft_plan &my_plan, fft_complex &fft_of_rhogrid)
{
  double k_interval = 1.005;
  int num_kbins     = (int)(log(sqrt(3.) * PMGRID / 0.95) / log(k_interval));
  if(All.NeutrinoScheme == 4.0)
    {
      int b;
      double start_k = 2. * MPI * 0.95 / (All.BoxSize / 1e3);
      double kk;
      double *pk_b, *pk_nub;
      pk_b   = (double *)malloc(num_kbins * sizeof(double));
      pk_nub = (double *)malloc(num_kbins * All.NNeutrino * sizeof(double));

      /* initialise b arrays */
      k_array0[0] = start_k;
      for(b = 1; b < (num_kbins + 1); b++)
        {
          k_array0[b] = k_array0[b - 1] * k_interval;
        }

      for(b = 0; b < num_kbins; b++)
        {
          pk_b[b] = 0.;
          for(i = 0; i < All.NNeutrino; i++)
            {
              pk_nub[b * All.NNeutrino] = 0.;
            }
          count_b[b] = 0.;
          k_array[b] = sqrt(k_array0[b] * k_array0[b + 1]);
        }

      for(i = 0; i < output_time_size; i++)
        {
          if(All.Time >= output_time_array[i] && All.A_Last_PM_Step < output_time_array[i + 1])
            {
              if(ThisTask == 0)
                {
                  printf("here the time is an output\n");
                  int nj;
                  FILE *fp;
                  char nu_txt[300];

                  sprintf(nu_txt, "%s_%d.txt", All.Nu_Pk_Txt, i);
                  fp = fopen(nu_txt, "w");

                  double nu_temp;
                  for(nj = 0; nj < num_kbins; nj++)
                    {
                      if(old_pk_b[nj] > 1e-7)
                        {
                          nu_temp = 0.;
                          for(j = 0; j < All.NNeutrino; j++)
                            {
                              nu_temp += old_pk_nu_b[nj * All.NNeutrino + j] * old_pk_nu_b[nj * All.NNeutrino + j];
                            }
                          nu_temp /= (double)All.NNeutrino;
                          fprintf(fp, "%f\t %.20f\n", k_array[nj], nu_temp / (old_pk_b[nj] * old_pk_b[nj]));
                        }
                    }
                  fclose(fp);
                  printf("finished printing nu_pk.txt no.%d time %f\n", i, output_time_array[i]);
                }
            }
        }

      /* for the 4.0 version correction, we need to calculate the cdm pk first, so need to loop twice: count number
       * in each bin; calculate pk;
       */
      for(x = 0; x < GRIDX; x++)
        for(y = myplan.slabstart_y; y < myplan.slabstart_y + myplan.nslab_y; y++)
          for(z = 0; z < GRIDz; z++)
            {
              int xx, yy, zz;

              if(x >= (GRIDX / 2))
                xx = x - GRIDX;
              else
                xx = x;
              if(y >= (GRIDY / 2))
                yy = y - GRIDY;
              else
                yy = y;
              if(z >= (GRIDZ / 2))
                zz = z - GRIDZ;
              else
                zz = z;

              double kx = kfacx * xx;
              double ky = kfacy * yy;
              double kz = kfacz * zz;

              double k2 = kx * kx + ky * ky + kz * kz;
              kk        = pow(k2, 0.5) * 2. * M_PI * 1e3 / All.Boxsize;

#ifndef FFT_COLUMN_BASED
              large_array_offset ip = ((large_array_offset)GRIDz) * (GRIDX * (y - myplan.slabstart_y) + x) + z;
#endif
              for(b = 0; b < num_kbins; b++)
                {
                  if(kk >= k_array0[b] && kk < k_array0[b + 1])
                    {
                      count_b[b] = count_b[b] + 1.;
                      pk_b[b] += (fft_of_rhogrid[ip].re * fft_of_rhogrid[ip].re + fft_of_rhogrid[ip].im * fft_of_rhogrid[ip].im);
                    }
                }
            }
      if(All.NumCurrentTiStep == 0)
        {
          int rd;
          double ratio_temp = 1.;
          for(b = 0; b < num_kbins; b++)
            {
              if(count_b[b] > 0.)
                {
                  old_pk_b[b] = sqrt(pk_b[b] / count_b[b]);

                  for(b = 0; b < rd_size_cal; b++)
                    {
                      if(k_array[b] >= rd_array_k[rd] && k_array[b] < rd_array_k[rd + 1])
                        {
                          ratio_temp = (rd_array_pk[rd] + rd_array_pk[rd + 1]) / 2.;
                        }
                    }
                  for(i = 0; i < All.NNeutrino; i++)
                    {
                      old_pk_nu_b[b * All.NNeutrino + i] = old_pk_b[b] * sqrt[ratio_temp];
                    }
                }
            }
          All.A_Last_PM_Step = All.Time;
        }
      double fnu_total, roneu_temp_total;
      double *fnu, *roneu_temp;
      fnu              = (double *)malloc(All.NNeutrino * sizeof(double));
      roneu_temp       = (double *)malloc(All.NNeutrino * sizeof(double));
      roneu_temp_total = 0.;
      for(i = 0; i < All.NNeutrino; i++)
        {
          roneu_temp[i] = neutrino_integration(All.Time, All.mass[i], All.xi[i]);
          roneu_temp_total += roneu_temp[i];
        }
      fnu_total = 0.;
      for(i = 0; i < All.NNeutrino; i++)
        {
          fnu[i] = roneu_temp[i] / (roneu_temp_total + (All.Omega0 - All.Omega_nu0_frstr) / pow(All.Time, 3));
          fnu_total += fnu[i];
        }

      if(All.NumCurrentTiStep > 0)
        {
          double a_spacing;
          a_spacing = (All.Time - All.A_Last_PM_Step) / (double)(All.FrstrInterval);

          double *a_inte_series;
          double *s_inte_series;
          a_inte_series = (double *)malloc((All.FrstrInterval + 1) * sizeof(double));
          s_inte_series = (double *)malloc((All.FrstrInterval + 1) * sizeof(double));

          for(i = 0; i <= All.FrstrInterval; i++)
            {
              a_inte_series[i] = All.A_Last_PM_Step + i * a_spacing;
            }
          // time series calculation finished
          for(i = 1; i <= All.FrstrInterval; i++)
            {
              s_inte_series[i] =
                  a_to_s(i, a_inte_series, All.A_Last_PM_Step, All.Time) * c /
                  mpc_to_m;  // this unit is because in later frstr phi integration we removed this unit from the s series
              // printf("time %f i %d\n", s_inte_series[i], i);
            }
          s_inte_series[0] = 0.;

          for(b = 0; b < num_kbins; b++)
            {
              if(count_b[b] > 0)
                {
                  pk_b[b] = sqrt(pk_b[b] / count_b[b]);
                  for(i = 0; i < All.NNeutrino; i++)
                    {
                      pk_nub[b * All.NNeutrino + i] =
                          frstr(k_array[b], old_pk_nu_b[b * All.NNeutrino + i], old_pk_b[b], pk_b[b], a_inte_series, s_inte_series,
                                All.A_Last_PM_Step, All.Time, All.NuMass[i], All.Xi[i]);
                    }
                }
            }

          // printf("fnu %f xi %f--------------\n", fnu, All.xi_1);

          for(x = 0; x < GRIDX; x++)
            for(y = myplan.slabstart_y; y < myplan.slabstart_y + myplan.nslab_y; y++)
              for(z = 0; z < GRIDz; z++)
                {
                  int xx, yy, zz;

                  if(x >= (GRIDX / 2))
                    xx = x - GRIDX;
                  else
                    xx = x;
                  if(y >= (GRIDY / 2))
                    yy = y - GRIDY;
                  else
                    yy = y;
                  if(z >= (GRIDZ / 2))
                    zz = z - GRIDZ;
                  else
                    zz = z;

                  double kx = kfacx * xx;
                  double ky = kfacy * yy;
                  double kz = kfacz * zz;

                  double k2 = kx * kx + ky * ky + kz * kz;
                  kk        = pow(k2, 0.5) * 2. * M_PI * 1e3 / All.Boxsize;

#ifndef FFT_COLUMN_BASED
                  large_array_offset ip = ((large_array_offset)GRIDz) * (GRIDX * (y - myplan.slabstart_y) + x) + z;
#endif

                  int b0 = 0;
                  for(b = 0; b < num_kbins; b++)
                    {
                      if(kk >= k_array0[b] && kk < k_array0[b + 1])
                        {
                          double sum = 0.;
                          for(i = 0; i < All.NNeutrino; i++)
                            {
                              sum += fnu[i] * pk_nub[b * All.NNeutrino + i];
                            }
                          fft_of_rhogrid[ip].re =
                              fft_of_rhogrid[ip].re * (1. - fnu_total) + fft_of_rhogrid[ip].re * fabs(sum / pk_b[b]);
                          fft_of_rhogrid[ip].im =
                              fft_of_rhogrid[ip].im * (1. - fnu_total) + fft_of_rhogrid[ip].im * fabs(sum / pk_b[b]);
                          b0 = b;
                        }
                    }
                }

          for(b = 0; b < num_kbins; b++)
            {
              for(i = 0; i < All.NNeutrino; i++)
                {
                  old_pk_nu_b[b * All.NNeutrino + i] = pk_nub[b * All.NNeutrino + i];
                }
              old_pk_b[b] = pk_b[b];
            }

          All.A_Last_PM_Step = All.Time;
        }

      if(ThisTask == 0)
        {
          printf("time now %f time max %f\n", All.Time, All.TimeMax);
        }

      if(fabs(All.Time - All.TimeMax) < 1e-6)
        {
          if(ThisTask == 0)
            {
              printf("here the time is time max\n");
              int nj;
              FILE *fp;
              char nu_txt[300];
              double nu_temp;

              sprintf(nu_txt, "%s_%d.txt", All.Nu_Pk_Txt, (output_time_size));
              fp = fopen(nu_txt, "w");

              for(nj = 0; nj < num_kbins; nj++)
                {
                  if(old_pk_b[nj] > 1e-7)
                    {
                      nu_temp = 0.;
                      for(i = 0; i < All.NNeutrino; i++)
                        {
                          nu_temp += old_pk_nu_b[nj * All.NNeutrino + i] * old_pk_nu_b[nj * All.NNeutrino + i];
                        }
                      nu_temp /= (double)All.NNeutrino;

                      fprintf(fp, "%f\t %.20f\n", k_array[nj], nu_temp / (old_pk_b[nj] * old_pk_b[nj]));
                    }
                }
              fclose(fp);

              printf("finished printing nu_pk.txt\n");
            }
        }
    }
}
