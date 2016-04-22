#ifndef PROJECT2_REPULSIVEGAUSSIAN_H
#define PROJECT2_SIMPLEGAUSSIAN_H
#include "wavefunction.h"

class RepulsiveGaussian : public WaveFunction {
public:
    RepulsiveGaussian(class System* system, double alpha, double beta, double a);
    double evaluate(std::vector<class Particle*> particles);
    double computeDoubleDerivative(std::vector<class Particle*> particles);
    double computeDerivativeWrtAlpha(std::vector<Particle *> particles);
    double computeDerivativeWrtBeta(std::vector<Particle *> particles);
    std::vector<double> computeDerivative(std::vector<class Particle*> particles);

private:
    double m_a = 0;
};

#endif // PROJECT2_SIMPLEGAUSSIAN_H
