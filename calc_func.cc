#include <iostream>
#include "coulcorr_param.h"
using namespace std;

int main()
{
  double param[3] = {1,5,1.1}; //lambda, R, alpha
  for(double q=0; q<0.2; q+=0.002)
  {
    cout << q << "\t" << coulcorr_param(&q,&param[0]) << endl;
  }
}