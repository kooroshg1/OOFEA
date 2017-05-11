// local includes
#include "continuum_mechanics.h"

namespace Extra {
    ElasticityConstruction::ElasticityConstruction() {

    }

    ElasticityConstruction::~ElasticityConstruction() {

    }
    // Kronecker delta function
    inline double ElasticityConstruction::kronecker_delta(unsigned int i,
                                unsigned int j) {
      return i == j ? 1. : 0.;
    }

    double ElasticityConstruction::elasticity_tensor(double young_modulus,
                                                     double poisson_ratio,
                                                     unsigned int i,
                                                     unsigned int j,
                                                     unsigned int k,
                                                     unsigned int l) {
      // Define the Poisson ratio and Young's modulus
//      const double nu = 0.0;
//      const double E  = 1.;

      // Define the Lame constants (lambda_1 and lambda_2) based on nu and E
      const double lambda = young_modulus * poisson_ratio / ((1. + poisson_ratio) * (1. - 2.*poisson_ratio));
      const double mu = 0.5 * young_modulus / (1. + poisson_ratio);

      return lambda * this->kronecker_delta(i, j) * this->kronecker_delta(k, l) +
             mu * (this->kronecker_delta(i, k) * this->kronecker_delta(j, l) +
                   this->kronecker_delta(i, l) * this->kronecker_delta(j, k));
    }
}
