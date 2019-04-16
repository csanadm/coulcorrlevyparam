#include <cmath>
#define HBARC 0.197327
const double MASS_PI = 0.13957;
const double ALPHAEM = 0.007297;

double Gamow(const double *x, const double *par)
{
        double eta = ALPHAEM*MASS_PI/(x[0]); //x[0] = Q = 2k
        return 1./( M_PI*eta*2.0/(exp(2.0*M_PI*eta)-1) );
}

double coulcorr_param(const double *x, const double *par)
{
  double q = x[0];
  double R = par[1];
  double alpha = par[2];

  double aA =  0.213918;
  double aB = -0.503681;
  double aC =  0.0292535;
  double aD = -0.912421;
  double aE =  0.00183488;

  double bA =  1.54117;
  double bB =  0.38916;
  double bC =  1.95399;
  double bD = -1.35931;
  double bE = -0.061474;
  double bF =  0.770834;

  double cA = -1.23221;
  double cB =  0.00565967;
  double cC =  0.227756;
  double cD =  0.0086754;
  double cE =  0.00351454;
  double cF =  1.72461;

  double dA = -0.053639;
  double dB =  0.454527;
  double dC = -0.112904;
  double dD =  1.58131;
  double dE =  0.890556;
  double dF = -8.49402;

  double a_AR = pow((aA * alpha + aB),2) + pow((aC * R + aD),2) + aE * pow((alpha * R + 1),2);
  double b_AR = (1. + bA*pow(R,bB) -pow(alpha,bC)) / (alpha*alpha*R*(pow(alpha,bD)+bE*pow(R,bF)));
  double c_AR = (cA  + pow(alpha,cB) + cC*pow(R,cD))/cE*pow(alpha/R,cF);
  double d_AR = (dA + (pow(R,dB) + dC*pow(alpha,dF))/(pow(R,dD)*pow(alpha,dE)));
  double parametrization = 1./(1. + 0.03243*a_AR*R/pow(alpha,1.)/(1.+b_AR*10.1355*R/pow(alpha,1.)*q+c_AR*pow(10.1355*R/pow(alpha,1.)*q,2)+pow(d_AR*10.1355*R/pow(alpha,1.)*q,4)));

  double Aa =  0.169458;
  double Ab =  0.00908016;
  double Ac = -0.0193868;
  double Ad = -0.00747829;
  double Ae =  0.000952962;
  double Af =  0.000227839;

  double Ba = 21.8577;
  double Bb =  5.14268;
  double Bc =  0.955366;
  double Bd = -0.534315;
  double Be = -0.0339711;
  double Bf =  0.00299237;

  double A = Aa + Ab*alpha + Ac*R + Ad*alpha*R + Ae*R*R + Af*alpha*alpha*R*R;
  double B = Ba + Bb*alpha + Bc*R + Bd*alpha*R + Be*R*R + Bf*alpha*alpha*R*R;
  double exponential_smoothing = 1. + A * exp( - B * q );

  double q0 = 0.08;
  double Dq = 0.025;
  double cutoff = 1. / (1.+exp((q-q0)/Dq));

  return (1. - cutoff) * exponential_smoothing + cutoff * Gamow(x,par) * parametrization;
}
