#ifndef PROJECT2_HARMONICOSCILLATORELECTRONS_H
#define PROJECT2_HARMONICOSCILLATORELECTRONS_H
#include "hamiltonian.h"
#include <vector>

class HarmonicOscillatorElectrons : public Hamiltonian {
public:
    HarmonicOscillatorElectrons(System* system, double omega, bool analyticalKinetic);
    double computeLocalEnergy(std::vector<Particle*> particles);

private:
    double m_omega = 0;
};

#endif // PROJECT2_HARMONICOSCILLATORELECTRONS_H
