// This code reproduces the Coulomb correction based on numerical calculations with a Levy source
// Version: 2019/05/10
// Authors: S. Lokos and M. Csanad
#include <cmath>
#include "coulcorr_param.h"

using namespace std;

#define HBARC 0.1973269788       // hbar*c in GeV*fm 
#define PI_MASS 0.13957018       // Pion mass in GeV/c^2
#define ALPHAEM 0.0072973525664  // 1.0/137.035999139, fine structure constant
#define PREFACTOR 0.016215165    // ALPHAEM * PI * PI_MASS / HBARC
#define SQR(x) ((x)*(x))

// This returns the Coulomb correction, parametrized for a fully chaotic, single component (no-halo) source with Levy scale R and Levy index alpha
double coulcorr_param(const double *x, const double *par)
{
  double q = x[0];
  double R = par[1];
  double alpha = par[2];

  double aA =  0.360596;
  double aB = -0.545079;
  double aC =  0.0347488;
  double aD = -1.30389;
  double aE =  0.00377961;

  double bA =  2.04017;
  double bB =  0.559724;
  double bC =  2.47224;
  double bD = -1.26815;
  double bE = -0.117674;
  double bF =  0.527379;

  double cA = -1.00015;
  double cB =  0.000123187;
  double cC =  7.98451e-05;
  double cD =  0.269861;
  double cE =  2.56089e-05;
  double cF =  1.75202;

  double dA = 0.00263107;
  double dB =  -0.131239;
  double dC = -0.831485;
  double dD =  1.57528;
  double dE =  0.275676;
  double dF = 0.0493686;

  double a_AR = SQR((aA * alpha + aB)) + SQR((aC * R + aD)) + aE * SQR((alpha * R + 1));
  double b_AR = (1. + bA*pow(R,bB) -pow(alpha,bC)) / (alpha*alpha*R*(pow(alpha,bD)+bE*pow(R,bF)));
  double c_AR = (cA  + pow(alpha,cB) + cC*pow(R,cD))/cE*pow(alpha/R,cF);
  double d_AR = (dA + (pow(R,dB) + dC*pow(alpha,dF))/(pow(R,dD)*pow(alpha,dE)));
  double t = R*q/HBARC/alpha;
  double parametrization = 1./(1. + PREFACTOR*a_AR*R/alpha/(1.+b_AR*t+c_AR*t*t+d_AR*t*t*t*t));

  double Aa =  0.207369;
  double Ab = -0.0099871;
  double Ac = -0.0267108;
  double Ad = -0.00372829;
  double Ae =  0.00118981;
  double Af =  0.000159437;

  double Ba = 25.805;
  double Bb =  4.01674;
  double Bc =  0.0087332;
  double Bd = -0.256058;
  double Be =  0.0107724;
  double Bf = -0.00270474;

  double A = Aa + Ab*alpha + Ac*R + Ad*alpha*R + Ae*R*R + Af*alpha*alpha*R*R;
  double B = Ba + Bb*alpha + Bc*R + Bd*alpha*R + Be*R*R + Bf*alpha*alpha*R*R;
  double exponential_smoothing = 1. + A * exp( - B * q );

  double q0 = 0.08;
  double Dq = 0.02;
  double cutoff = 1. / (1.+exp((q-q0)/Dq));

  return 1. / ( (1. - cutoff) * exponential_smoothing + cutoff * Gamow(x,par) * parametrization );
}

// Inverse of the usual Gamow factor
double Gamow(const double *x, const double *par)
{
        double eta = ALPHAEM*PI_MASS/(x[0]); //x[0] is Q = 2k = |p1-p2|
        return 1./( M_PI*eta*2.0/(exp(2.0*M_PI*eta)-1) );
}
