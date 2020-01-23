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
   
    double aA =  0.26984;
    double aB = -0.491231;
    double aC =  0.0352313;
    double aD = -1.31628;
    double aE =  0.00359148;

    double bA =  2.37267;
    double bB =  0.586309;
    double bC =  2.24867;
    double bD = -1.43278;
    double bE = -0.0521642;
    double bF =  0.729434;

    double cA = -4.30347 ;
    double cB =  1.17849e-05;
    double cC =  3.30346;
    double cD =  1.27273e-06 ;
    double cE =  3.03399e-06;
    double cF =  1.68883;

    double dA =  0.000568486;
    double dB = -0.805271;
    double dC = -0.192606 ;
    double dD =  2.77504;
    double dE =  2.02951;
    double dF =  1.07906; 
   
    double a_AR = pow((aA * alpha + aB),2) + pow((aC * R + aD),2) + aE * pow((alpha * R + 1),2);
    double b_AR = (1. + bA*pow(R,bB) -pow(alpha,bC)) / (alpha*alpha*R*(pow(alpha,bD)+bE*pow(R,bF)));
    double c_AR = (cA  + pow(alpha,cB) + cC*pow(R,cD))/cE*pow(alpha/R,cF);
    double d_AR = (dA + (pow(R,dB) + dC*pow(alpha,dF))/(pow(R,dD)*pow(alpha,dE)));
    double t = R*q/HBARC/alpha;
    double parametrization = 1./(1. + PREFACTOR*a_AR*R/alpha/(1.+b_AR*t+c_AR*t*t+d_AR*t*t*t*t));
       
    double Aa =  0.126253;
    double Ab =  0.053848;
    double Ac = -0.00912627;
    double Ad = -0.018459;
    double Ae =  0.000851755;
    double Af =  0.000417179;

    double Ba = 19.3162;
    double Bb =  5.58961;
    double Bc =  2.26264;
    double Bd = -1.28486;
    double Be = -0.0821616;
    double Bf =  0.0238446; 
   
    double A = Aa + Ab*alpha + Ac*R + Ad*alpha*R + Ae*R*R + Af*alpha*alpha*R*R;
    double B = Ba + Bb*alpha + Bc*R + Bd*alpha*R + Be*R*R + Bf*alpha*alpha*R*R;
    double exponential_smoothing = 1. + A * exp( - B * q );
   
   
    double q0 = 0.07;
    double n = 20.0;
    double cutoff = 1. / ( 1.+pow(q/q0,n) );
   
   
    return 1./(  (1. - cutoff) * exponential_smoothing + cutoff * Gamow(x,par) * parametrization);
}

// Inverse of the usual Gamow factor
double Gamow(const double *x, const double *par)
{
    double eta = ALPHAEM*PI_MASS/(x[0]); //x[0] is Q = 2k = |p1-p2|
    return 1./( M_PI*eta*2.0/(exp(2.0*M_PI*eta)-1) );
}
