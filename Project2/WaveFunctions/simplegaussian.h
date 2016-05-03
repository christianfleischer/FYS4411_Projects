#ifndef PROJECT2_SIMPLEGAUSSIAN_H
#define PROJECT2_SIMPLEGAUSSIAN_H
#include "wavefunction.h"

class SimpleGaussian : public WaveFunction {
public:
    SimpleGaussian(class System* system, double alpha);
    double evaluate(std::vector<class Particle*> particles);
    double computeDoubleDerivative(std::vector<class Particle*> particles);
    double computeDerivativeWrtAlpha(std::vector<Particle *> particles);
    double computeDerivativeWrtBeta(std::vector<Particle *> particles);
    double computeMetropolisRatio(std::vector<Particle *> particles, int randomParticle,
                                  std::vector<double> positionChange);
    std::vector<double> computeDerivative(std::vector<class Particle*> particles);
    void updateSlaterDet(int randomParticle);
};

#endif // PROJECT2_SIMPLEGAUSSIAN_H
