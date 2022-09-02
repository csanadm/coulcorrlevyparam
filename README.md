# Coulomb correction parametrization based on a Levy source

## Description
This package contains a calculation for quantum-statistical correlation functions, including Coulomb-correction, interpolated from numerical integration calculations, detailed in the below mentioned calculations. Range of validity is R=2-12 fm for the HBT (Lévy) scale and alpha=0.8-2.0 for the Lévy exponent. The lambda parameter is included only in the correlation function calculations, as done in a Bowler-Sinyukov type of treatment.

## File content
- [**README.md**](https://github.com/csanadm/coulcorrlevyparam/blob/master/README.md): This README file
- [**Makefile**](https://github.com/csanadm/coulcorrlevyparam/blob/main/Makefile): Using `make all`, it will create an executable
- [**coulcorr_param.h**](https://github.com/csanadm/coulcorrlevyparam/blob/main/coulcorr_param.h): This contains the function declarations for the calculation
- [**coulcorr_param.cc**](https://github.com/csanadm/coulcorrlevyparam/blob/main/coulcorr_param.cc): This contains the calculations themselves
- [**calc_func.cc**](https://github.com/csanadm/coulcorrlevyparam/blob/main/calc_func.cc): Example calculation of the correlation function, including Coulomb correction

## Publications
- Expanded empirical formula for Coulomb final state interaction in the presence of Lévy sources
Máté Csanád, Sándor Lökös, Márton Nagy
Phys.Part.Nucl. 51 (2020) 3, 238-242, [arXiv:1910.02231 [hep-ph]](https://arxiv.org/abs/1910.02231).
- Coulomb final state interaction in heavy ion collisions for Lévy sources
Máté Csanád, Sándor Lökös, Márton Nagy
Universe 5 (2019) 6, 133, [arXiv:1905.09714 [nucl-th]](https://arxiv.org/abs/1905.09714).
