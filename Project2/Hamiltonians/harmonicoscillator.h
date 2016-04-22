#ifndef PROJECT2_HARMONICOSCILLATOR_H
#define PROJECT2_HARMONICOSCILLATOR_H
#include "hamiltonian.h"
#include <vector>

class HarmonicOscillator : public Hamiltonian {
public:
    HarmonicOscillator(System* system, double omega, bool analyticalKinetic);
    double computeLocalEnergy(std::vector<Particle*> particles);

private:
    double m_omega = 0;
};

#endif // PROJECT2_HARMONICOSCILLATOR_H
