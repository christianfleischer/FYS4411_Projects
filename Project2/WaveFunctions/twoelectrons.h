#pragma once
#include "wavefunction.h"

class TwoElectrons : public WaveFunction {
public:
    TwoElectrons(class System* system, double alpha, double beta, double omega, double a, double C);
    double evaluate(std::vector<class Particle*> particles);
    double computeDoubleDerivative(std::vector<class Particle*> particles);
    double computeDerivativeWrtAlpha(std::vector<Particle *> particles);
    double computeDerivativeWrtBeta(std::vector<Particle *> particles);
    std::vector<double> computeDerivative(std::vector<class Particle*> particles);

private:
    double m_omega = 0;
    double m_a = 0;
    double m_C = 0;
};

