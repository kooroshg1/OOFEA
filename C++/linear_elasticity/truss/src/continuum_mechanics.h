#ifndef CONTINUUM_MECHANICS_H
#define CONTINUUM_MECHANICS_H

namespace Extra {
    class ElasticityConstruction {
    public:
        // Kronecker delta function
        ElasticityConstruction();
        ~ElasticityConstruction();

        // Rank-4 tensor for elasticity
        double elasticity_tensor(double young_modulus,
                                 double poisson_ratio,
                                 unsigned int i,
                                 unsigned int j,
                                 unsigned int k,
                                 unsigned int l);
    
        inline double kronecker_delta(unsigned int i,
                                      unsigned int j);
    };
}
#endif
