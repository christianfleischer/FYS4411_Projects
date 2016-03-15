#pragma once
#include "hamiltonian.h"
#include <vector>

class HarmonicOscillatorRepulsive : public Hamiltonian {
public:
    HarmonicOscillatorRepulsive(System* system, double omega, double a, double gamma, bool analyticalKinetic);
    double computeLocalEnergy(std::vector<Particle*> particles);

private:
    double m_omega = 0;
    double m_a = 0;
    double m_gamma = 0;
};
