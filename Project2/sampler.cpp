#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "sampler.h"
#include "system.h"
#include "particle.h"
#include "Hamiltonians/hamiltonian.h"
#include "WaveFunctions/wavefunction.h"

using namespace std;

Sampler::Sampler(System* system) {
    m_system = system;
    m_stepNumber = 0;
}

void Sampler::setNumberOfMetropolisSteps(int steps) {
    m_numberOfMetropolisSteps = steps;
}

void Sampler::setEnergy(double energy) {
    m_energy = energy;
}

void Sampler::setVariance(double variance) {
    m_variance = variance;
}

void Sampler::setAcceptanceRate(double acceptanceRate) {
    m_acceptanceRate = acceptanceRate;
}

void Sampler::setMeanDistance(double meanDistance) {
    m_meanDistance = meanDistance;
}

void Sampler::sample(bool acceptedStep, bool saveEnergies, bool savePositions) {
    // Make sure the sampling variable(s) are initialized at the first step.
    if (m_stepNumber == 0) {
        m_cumulativeEnergy =             0;
        m_cumulativeSquaredEnergy =      0;
        m_cumulativeWFuncDerivativeAlpha =    0;
        m_cumulativeWFuncEnergyAlpha =        0;
        m_cumulativeAcceptedSteps =      0;
    }

    /* Here you should sample all the interesting things you want to measure.
     * Note that there are (way) more than the single one here currently.
     */

    // Sample local energy
    double localEnergy = m_system->getHamiltonian()->
                         computeLocalEnergy(m_system->getParticles());
    m_cumulativeEnergy  += localEnergy;
    m_cumulativeSquaredEnergy += localEnergy*localEnergy;

    std::vector<Particle*> particles = m_system->getParticles();

    // Sample mean distance
    double r12 = 0;
    int numberOfParticles = m_system->getNumberOfParticles();
    int numberOfDimensions = m_system->getNumberOfDimensions();
    if (numberOfParticles == 2){
        std::vector<double> r1 = particles[0]->getPosition();
        std::vector<double> r2 = particles[1]->getPosition();
        for (int i=1; i < numberOfDimensions; i++){
            r12 += (r1[i]-r2[i])*(r1[i]-r2[i]);
        }
    }
    r12 = sqrt(r12);
    m_cumulativeDistance += r12;

    // Sample things needed for the steepest descent method:
    double waveFunction = m_system->getWaveFunction()->evaluate(particles);

    double waveFuncDerivativeAlpha = m_system->getWaveFunction()->computeDerivativeWrtAlpha(particles);
    waveFuncDerivativeAlpha /= waveFunction;
    m_cumulativeWFuncDerivativeAlpha += waveFuncDerivativeAlpha;
    m_cumulativeWFuncEnergyAlpha += waveFuncDerivativeAlpha*localEnergy;

    double waveFuncDerivativeBeta = m_system->getWaveFunction()->computeDerivativeWrtBeta(particles);
    waveFuncDerivativeBeta /= waveFunction;
    m_cumulativeWFuncDerivativeBeta += waveFuncDerivativeBeta;
    m_cumulativeWFuncEnergyBeta += waveFuncDerivativeBeta*localEnergy;

    // Sample whether the step is accepted or not in order to find total acceptance ratio
    m_cumulativeAcceptedSteps += acceptedStep;

    saveToFile(localEnergy, saveEnergies, savePositions);
    m_stepNumber++;
}

void Sampler::printOutputToTerminal() {
    int     np = m_system->getNumberOfParticles();
    int     nd = m_system->getNumberOfDimensions();
    int     ms = m_system->getNumberOfMetropolisSteps();
    int     p  = m_system->getWaveFunction()->getNumberOfParameters();
    double  ef = m_system->getEquilibrationFraction();
    double  ct = m_system->getComputationTime();
    std::vector<double> pa = m_system->getWaveFunction()->getParameters();

    cout << "    " << endl;
    cout << "  -- System info -- " << endl;
    cout << " Number of particles  : " << np << endl;
    cout << " Number of dimensions : " << nd << endl;
    cout << " Number of Metropolis steps run : 10^" << std::log10(ms) << endl;
    cout << " Number of equilibration steps  : 10^" << std::log10(std::round(ms*ef)) << endl;
    cout << endl;
    cout << "  -- Wave function parameters -- " << endl;
    cout << " Number of parameters : " << p << endl;
    for (int i=0; i < p; i++) {
        cout << " Parameter " << i+1 << " : " << pa.at(i) << endl;
    }
    cout << endl;
    cout << "  -- Results -- " << endl;
    cout << " Energy : " << m_energy << endl;
    cout << " Variance : " << m_variance << endl;
    cout << " Acceptance Rate : " << m_acceptanceRate << endl;
    if (m_meanDistance != 0) cout << " Mean Distance : " << m_meanDistance << endl;
    cout << endl;
    cout << "Computation Time : " << ct << endl;
    cout << endl;
}

void Sampler::computeAverages() {
    /* Compute the averages of the sampled quantities. You need to think
     * thoroughly through what is written here currently; is this correct?  It should be now.
     */
    // We only sample after the equilibration steps so we need to average by the number of steps after equilibration.
    m_energy = m_cumulativeEnergy / (double) m_stepNumber;
    m_squaredEnergy = m_cumulativeSquaredEnergy / (double) m_stepNumber;
    m_variance = m_squaredEnergy - m_energy*m_energy;

    m_waveFuncDerivativeAlpha = m_cumulativeWFuncDerivativeAlpha / (double) m_stepNumber;
    m_waveFuncEnergyAlpha = m_cumulativeWFuncEnergyAlpha / (double) m_stepNumber;

    m_waveFuncDerivativeBeta = m_cumulativeWFuncDerivativeBeta / (double) m_stepNumber;
    m_waveFuncEnergyBeta = m_cumulativeWFuncEnergyBeta / (double) m_stepNumber;

    m_acceptanceRate = m_cumulativeAcceptedSteps / (double) m_stepNumber;

    m_meanDistance = m_cumulativeDistance / (double) m_stepNumber;
}

void Sampler::saveToFile(double localEnergy, bool saveEnergies, bool savePositions){
    if (saveEnergies){
        if (m_stepNumber == 0) system("rm energies.dat");

        FILE *outfile = fopen("energies.dat", "a");
        fprintf(outfile, "%f\n", localEnergy);
        fclose(outfile);
    }

    if (savePositions){
        if (m_stepNumber == 0) system("rm positions.dat");

        FILE *outfile = fopen("positions.dat", "a");
        std::vector<Particle*> particles = m_system->getParticles();
        int numberOfParticles = m_system->getNumberOfParticles();
        int numberOfDimensions = m_system->getNumberOfDimensions();

        for (int i=0; i < numberOfParticles; i++){
            std::vector<double> r_i = particles[i]->getPosition();
            for (int j=0; j < numberOfDimensions; j++){
                fprintf(outfile, "%f ", r_i[j]);
            }
            fprintf(outfile, "\n");
        }
        fclose(outfile);
    }
}




