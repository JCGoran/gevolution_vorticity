#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>

/**
    differential equation for the growth rate D_1
**/

static double g1(double a, void *p)
{
  struct cosmology par=*(struct cosmology *) p;
  double g1 = 1./(pow(a, 3)*pow(1 - par.Omega_m + par.Omega_m/pow(a, 3),1.5));
  return g1;
}

/**
    computes and stores all the background functions
**/

double D1_prime(
    struct cosmology par, double a
)
{
    gsl_function G1;
    G1.function = &g1, G1.params = &par;
    const double delta_a = 0.005;
    double a_minus = a*(1 - delta_a);
    double a_plus  = a*(1 + delta_a);
    gsl_integration_workspace *w;
    double prec = 1E-5;
    double result, error;
    w=gsl_integration_workspace_alloc(10000);
    gsl_integration_qag (&G1, 0, a_plus, 0, prec, 10000,
		         GSL_INTEG_GAUSS61, w, &result, &error);

    double result_plus = (5./2*par.Omega_m*sqrt(1 - par.Omega_m + par.Omega_m/pow(a_plus, 3)))*result;

    gsl_integration_qag (&G1, 0, a_minus, 0, prec, 10000,
                         GSL_INTEG_GAUSS61, w, &result, &error);

    double result_minus = (5./2*par.Omega_m*sqrt(1 - par.Omega_m + par.Omega_m/pow(a_minus, 3)))*result;

    double result_der = (result_plus -  result_minus)/(0.01*a);

    result = result_der*Hconf(a, 1, par)*a;

    gsl_integration_workspace_free(w);

    return result;
}
