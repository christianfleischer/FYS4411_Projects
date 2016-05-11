#ifndef PROJECT2_TWOELECTRONS_H
#define PROJECT2_TWOELECTRONS_H
#include "wavefunction.h"

class TwoElectrons : public WaveFunction {
public:
    TwoElectrons(class System* system, double alpha, double beta, double omega, double a, double C);
    double evaluate(std::vector<class Particle*> particles);
    double computeDoubleDerivative(std::vector<class Particle*> particles);
    double computeDerivativeWrtAlpha(std::vector<Particle *> particles);
    double computeDerivativeWrtBeta(std::vector<Particle *> particles);
    double computeMetropolisRatio(std::vector<Particle *> particles, int randomParticle,
                                  std::vector<double> positionChange);
    std::vector<double> computeDerivative(std::vector<class Particle*> particles,
                                          int randomParticle);
    void updateSlaterDet(int randomParticle);

private:
    double m_omega = 0;
    double m_a = 0;
    double m_C = 0;
};

#endif // PROJECT2_TWOELECTRONS_H
