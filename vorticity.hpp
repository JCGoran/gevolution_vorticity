/////////////////////////
// vorticity.hpp
//////////////////////////
//
// Function use to extract the vorticity field from gevolution
//
// Last modified: January 2019
//
//////////////////////////

#ifndef VORTICITY_HEADER
#define VORTICITY_HEADER

#include "prng_engine.hpp"
#include "d1_prime.hpp"

#include <gsl/gsl_spline.h>

using namespace std;
using namespace LATfield2;

#ifndef MAX_LINESIZE
#define MAX_LINESIZE 2048
#endif

#ifndef Cplx
#define Cplx Imag
#endif

using namespace std;
using namespace LATfield2;

// should be larger than maximum Ngrid
#ifndef HUGE_SKIP
#define HUGE_SKIP   65536
#endif

//////////////////////////
// Description:
//   Compute the rotational part of the velocity field, in Fourier space
//
//
// Arguments:
//
//   vRFT          reference to the Fourier image of the rotational part of the velocity
//   viFT          reference to the input Fourier image of the velocity field
//
// Returns:
//
//////////////////////////

void projectFTvelocityVR(Field<Cplx> & vRFT, Field<Cplx> & viFT)
{

  const int linesize = vRFT.lattice().size(1);
  int i;
  Real * gridk2;
  Cplx * kshift;
  Real * gridk;
  rKSite k(vRFT.lattice());
  Real k2;
  Cplx tmp(0., 0.);

  gridk2 = (Real *) malloc(linesize * sizeof(Real));
  kshift = (Cplx *) malloc(linesize * sizeof(Cplx));
  gridk = (Real *) malloc(linesize * sizeof(Real));

  for (i = 0; i < linesize; i++)
    {
      gridk[i] = (Real) linesize * sin(M_PI * 2.0 * (Real) i / (Real) linesize);
      kshift[i] = gridk2[i] * Cplx(cos(M_PI * (Real) i / (Real) linesize), -sin(M_PI * (Real) i / (Real) linesize));
      gridk2[i] = gridk[i]*gridk[i];
    }


  k.first();
  if (k.coord(0) == 0 && k.coord(1) == 0 && k.coord(2) == 0)
    {
      vRFT(k, 0) = Cplx(0.,0.);
      vRFT(k, 1) = Cplx(0.,0.);
      vRFT(k, 2) = Cplx(0.,0.);
      k.next();
    }
  for (; k.test(); k.next())
    {

      k2 = gridk2[k.coord(0)] + gridk2[k.coord(1)] + gridk2[k.coord(2)];
      if ((k.coord(0) == 0 || k.coord(0) == linesize/2)&&
          (k.coord(1) == 0 || k.coord(1) == linesize/2)&&
      (k.coord(2) == 0 || k.coord(2) == linesize/2)  ) {
    vRFT(k, 0) = Cplx(0.,0.);
    vRFT(k, 1) = Cplx(0.,0.);
    vRFT(k, 2) = Cplx(0.,0.);
      }
      else {
      tmp = (gridk[k.coord(0)] * viFT(k, 0) + gridk[k.coord(1)] * viFT(k, 1) + gridk[k.coord(2)] * viFT(k, 2)) / k2;

      vRFT(k, 0) = (viFT(k, 0) - gridk[k.coord(0)] * tmp);
      vRFT(k, 1) = (viFT(k, 1) - gridk[k.coord(1)] * tmp);
      vRFT(k, 2) = (viFT(k, 2) - gridk[k.coord(2)] * tmp);}
    }


  free(gridk2);
  free(kshift);

}

//////////////////////////
// projectFTvelocityTh
//////////////////////////
// Description:
//   Compute the diverge of the velocity in Fourier space
//
// Arguments:
//   thFT       reference to the Fourier image of the divergence of the velocity field
//   viFT       reference to the Fourier image of the velocity field
//
// Returns:
//
//////////////////////////

void projectFTvelocityTh(Field<Cplx> & thFT, Field<Cplx> & viFT)
{
  const int linesize = thFT.lattice().size(1);
  int i;
  Real * gridk2;
  Real * gridk;
  Cplx * kshift;
  rKSite k(thFT.lattice());
  Cplx tmp(0., 0.);

  gridk2 = (Real *) malloc(linesize * sizeof(Real));
  kshift = (Cplx *) malloc(linesize * sizeof(Cplx));
  gridk = (Real *) malloc(linesize * sizeof(Real));

  for (i = 0; i < linesize; i++)
    {
      gridk[i] = (Real) linesize * sin(M_PI * 2.0 * (Real) i / (Real) linesize);
      kshift[i] = gridk2[i] * Cplx(cos(M_PI * (Real) i / (Real) linesize), -sin(M_PI * (Real) i / (Real) linesize));
      gridk2[i] = gridk[i]*gridk[i];
    }


  for (k.first(); k.test(); k.next())
    {
      thFT(k) = Cplx(0.,1.)*(gridk[k.coord(0)] * viFT(k, 0) +
                 gridk[k.coord(1)] * viFT(k, 1) +
                 gridk[k.coord(2)] * viFT(k, 2) );
    }

  free(gridk2);
  free(kshift);
  free(gridk);
}

//////////////////////////
// compute_vi_past_rescaled
//////////////////////////
// Description:
//   Compute the velocity field as v^i = T^i_0/T^0_0, if a = 1 then vi = a v^i
//   If T^0_0 = 0 the velocity field is set to be the one at the previous time step,
//   rescaled as v^i(a) = v^i(a_past) a*Hconf(a) dD1/da (velocity method = rescaled past)
//
// Arguments:
//   viFT       reference to the velocity field
//   source     reference to the field source (a^3 T^0_0)
//   Bi         reference to the field Bi (a^4 T^0_i)
//   phi        reference to the field phi
//   chi        reference to the field chi
//   vi_past    reference to the velocity field at the previous time step
// Returns:
//
//////////////////////////

void compute_vi_past_rescaled(cosmology & cosmo, Field<Real> * vi, Field<Real> * source = NULL, double a = 1., double a_past = 1., Field<Real> * Ti0 = NULL, Field<Real> * vi_past = NULL)
{

  Site xvi(vi->lattice());

  Real rescale = D1_prime(cosmo, a)/D1_prime(cosmo, a_past)*a/a_past;

  for(xvi.first(); xvi.test(); xvi.next())
    {

      if ( (*source)(xvi) < 1.E-300) {(*vi)(xvi,0)= (*vi_past)(xvi,0)*rescale;}
      else {(*vi)(xvi,0) = (*Ti0)(xvi,0)/(*source)(xvi);}

      if ( (*source)(xvi) < 1.E-300) {(*vi)(xvi,1)= (*vi_past)(xvi,1)*rescale;}
      else {(*vi)(xvi,1) = (*Ti0)(xvi,1)/(*source)(xvi);}

      if ( (*source)(xvi) < 1.E-300) {(*vi)(xvi,2)= (*vi_past)(xvi,2)*rescale;}
      else {(*vi)(xvi,2) = (*Ti0)(xvi,2)/(*source)(xvi);}
    }
}

// Store the velocity field at each time step in vi_past
void store_vi(Field<Real> * vi_past, Field<Real> * vi = NULL)
{
  Site x(vi_past->lattice());

  if (vi != NULL)
  {
   for(x.first(); x.test(); x.next())
    {
      (*vi_past)(x,0) = (*vi)(x,0);
      (*vi_past)(x,1) = (*vi)(x,1);
      (*vi_past)(x,2) = (*vi)(x,2);
    }
  }
}

// Computes the norm of v_rotational (only used for snapshots!)
void compute_norm2_vR(
    	      Field<Real> * vR,
    	      Field<Real> * norm_vR2
    	      )
{
  Site x(norm_vR2->lattice());
  for(x.first(); x.test(); x.next()){
    (*norm_vR2)(x) = pow((*vR)(x, 0),2)+ pow((*vR)(x, 1),2)+ pow((*vR)(x, 2),2);
  }
}


// Computes the norm of the above (only used for snapshots!)
void compute_norm_w(
                    Field<Real> * norm2_vR,
                    Field<Real> * norm_w
                    )
{
  Site x(norm_w->lattice());
  for(x.first(); x.test(); x.next()){
    (*norm_w)(x) = 8.0*(*norm2_vR)(x)-((*norm2_vR)(x-0-1-2) + (*norm2_vR)(x-0-1+2) + (*norm2_vR)(x-0+1-2) + (*norm2_vR)(x+0-1-2) + \
    			       (*norm2_vR)(x-0+1+2) + (*norm2_vR)(x+0-1+2) + (*norm2_vR)(x+0+1-2) + (*norm2_vR)(x+0+1+2) );
    (*norm_w)(x) = sqrt(abs( (*norm_w)(x) ));

  }
}


//////////////////////////
// compute_sigma2_rescaled
//////////////////////////
// Description:
//   Compute the trace of the velocity dispertion tensor as sigma2 = -(T^1_1 + T^2_2 + T^3_3) / T^0_0 - <v>^2
//   If T^0_0 = 0 the velocity dispertion is set to be the one at the previous time step rescaled by linear velocity growth
//
// Arguments:
//   sigma2         reference to the velocity dispertion scalare
//   source         reference to the field source (a^3 T^0_0)
//   Sij            reference to the field Sij (a^3 T^i_j for the diagonal)
//   vi             reference to the velocity field
//   sigma2_past    reference to the velocity dispertion at the previous time step
// Returns:
//
//////////////////////////

void compute_sigma2_rescaled(cosmology & cosmo, Field<Real> * sigma2, Field<Real> * source = NULL, Field<Real> * Sij = NULL, Field<Real> * vi = NULL, Field<Real> * sigma2_past = NULL, double a = 1., double a_past = 1.)
{

  Site xsigma(sigma2->lattice());

  Real rescale = D1_prime(cosmo, a)/D1_prime(cosmo, a_past)*a/a_past;

  for(xsigma.first(); xsigma.test(); xsigma.next())
    {

      if ( (*source)(xsigma) < 1.E-300) {(*sigma2)(xsigma)= (*sigma2_past)(xsigma)*rescale*rescale;}
      else {(*sigma2)(xsigma) = ((*Sij)(xsigma, 0, 0) + (*Sij)(xsigma, 1, 1) + (*Sij)(xsigma, 2, 2) )/(*source)(xsigma)
      -((*vi)(xsigma,0)*(*vi)(xsigma,0) + (*vi)(xsigma,1)*(*vi)(xsigma,1) + (*vi)(xsigma,2)*(*vi)(xsigma,2));
           }

    }
}


// Store the velocity field at each time step in vi_past
void store_sigma2(Field<Real> * sigma2_past, Field<Real> * sigma2 = NULL)
{
  Site x(sigma2_past->lattice());

  if (sigma2 != NULL)
    {
      for(x.first(); x.test(); x.next())
    {
      (*sigma2_past)(x,0) = (*sigma2)(x,0);
      (*sigma2_past)(x,1) = (*sigma2)(x,1);
      (*sigma2_past)(x,2) = (*sigma2)(x,2);
    }
    }
}

#endif
