#include "hamiltonian.h"
#include "../system.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"

Hamiltonian::Hamiltonian(System* system, bool analyticalKinetic) {
    m_system = system;
    m_analyticalKinetic = analyticalKinetic;
}

double Hamiltonian::computeKineticEnergy(std::vector<Particle*> particles){
    // Compute the kinetic energy using numerical differentiation.

    double numberOfParticles = m_system->getNumberOfParticles();
    double numberOfDimensions = m_system->getNumberOfDimensions();
    double h = 1e-4;

    // Evaluate wave function at current step
    double waveFunctionCurrent = m_system->getWaveFunction()->evaluate(particles);
    double kineticEnergy = 0;

    for (int i=0; i < numberOfParticles; i++){
        for (int j=0; j < numberOfDimensions; j++){

            // Evaluate wave function at forward step
            particles[i]->adjustPosition(h, j);
            double waveFunctionPlus = m_system->getWaveFunction()->evaluate(particles);

            // Evaluate wave function at backward step
            particles[i]->adjustPosition(-2*h, j);
            double waveFunctionMinus = m_system->getWaveFunction()->evaluate(particles);

            // Part of numerical diff
            kineticEnergy -= (waveFunctionPlus - 2*waveFunctionCurrent + waveFunctionMinus);

            // Move particles back to original position
            particles[i]->adjustPosition(h, j);
        }
    }
    // Other part of numerical diff. Also divide by evaluation of current wave function
    // and multiply by 0.5 to get the actual kinetic energy.
    kineticEnergy = 0.5*kineticEnergy / (waveFunctionCurrent*h*h);
    return kineticEnergy;
}
