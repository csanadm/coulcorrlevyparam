#include <cstdlib>
#include <iostream>
#include "coulcorr_param.h"
using namespace std;

int main(int argc, char** argv)
{
  if(argc!=4)
  {
    cerr << "Usage: " << argv[0] << " <lambda> <R> <alpha>" << endl;
    return 1;
  }
  double param[3] = {atof(argv[1]),atof(argv[2]),atof(argv[3])}; //lambda, R, alpha
  for(double q=0.001; q<0.2; q+=0.001)
  {
    double lambda = param[0];
    double R = param[1];
    double alpha = param[2];
    double coulcorr = coulcorr_param(&q,&param[0]);
    double c2 = 1+exp(-pow(q*R/HBARC,alpha));
    cout << q << "\t" << 1-lambda+lambda*coulcorr*c2 << "\t" << coulcorr << endl;
  }
}
