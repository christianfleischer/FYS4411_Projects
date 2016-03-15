#include "system.h"
#include <cassert>
#include "sampler.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "Hamiltonians/hamiltonian.h"
#include "InitialStates/initialstate.h"
#include "Math/random.h"
#include <iostream>
#include <cmath>
#include <time.h>

using namespace std;

bool System::metropolisStep() {
    /* Perform the actual Metropolis step: Choose a particle at random and
     * change it's position by a random amount, and check if the step is
     * accepted by the Metropolis test (compare the wave function evaluated
     * at this new position with the one at the old position).
     */

    // Evaluate the wave function for current positions
    double waveFunctionOld = m_waveFunction->evaluate(m_particles);

    // Choose a random particle and change its position by a random amount creating a trial state
    int randomParticle = Random::nextInt(m_numberOfParticles);
    std::vector<double> positionChange(m_numberOfDimensions);

    for (int i=0; i<m_numberOfDimensions; i++){
        positionChange[i] = (Random::nextDouble()*2-1)*m_stepLength;
        m_particles[randomParticle]->adjustPosition(positionChange[i], i);
    }

    // Evaluate the wave function for the trial state
    double waveFunctionNew = m_waveFunction->evaluate(m_particles);
    // Metropolis ratio
    double qratio = waveFunctionNew*waveFunctionNew / (waveFunctionOld*waveFunctionOld);

    // Check if trial state is accepted
    if (Random::nextDouble() <= qratio){
        return true;
    }

    for (int i=0; i<m_numberOfDimensions; i++){
        // If trial state is not accepted, revert to old position for chosen particle (revert to old state)
        m_particles[randomParticle]->adjustPosition(-positionChange[i], i);
    }

    return false;
}

bool System::metropolisStepImpSampling(){

    // Choose a random particle to change the position of
    int randomParticle = Random::nextInt(m_numberOfParticles);
    std::vector<double> positionChange(m_numberOfDimensions);
    double D = 0.5;

    // Evaluate wave function for current state
    double waveFunctionOld = m_waveFunction->evaluate(m_particles);
    // Keep old position for Greens function
    std::vector<double> positionOld = m_particles[randomParticle]->getPosition();

    // Change position of random particle
    for (int i=0; i < m_numberOfDimensions; i++){
        positionChange[i] = Random::nextGaussian(0., sqrt(m_dt)); + D*m_dt*quantumForce(randomParticle)[i];
        m_particles[randomParticle]->adjustPosition(positionChange[i], i);
    }

    // Evaluate wave function for trial state
    double waveFunctionNew = m_waveFunction->evaluate(m_particles);
    // Keep new position for Greens function
    std::vector<double> positionNew = m_particles[randomParticle]->getPosition();

    // Evaluate Greens functions and find Metropolis-Hastings ratio:
    double GreensFunctionNew = calculateGreensFunction(randomParticle, positionNew, positionOld);

    for (int i=0; i < m_numberOfDimensions; i++){
        m_particles[randomParticle]->adjustPosition(-positionChange[i], i);
    }

    double GreensFunctionOld = calculateGreensFunction(randomParticle, positionOld, positionNew);

    double qratio = (GreensFunctionNew*waveFunctionNew*waveFunctionNew) / (GreensFunctionOld*waveFunctionOld*waveFunctionOld);

    // If move is accepted give the random particle the new position, otherwise keep the old position
    if (Random::nextDouble() <= qratio){
        for (int i=0; i<m_numberOfDimensions; i++){
            m_particles[randomParticle]->adjustPosition(positionChange[i], i);
        }
        return true;
    }

    return false;
}

std::vector<double> System::quantumForce(int particle){
    // Calculate the quantum (drift) force. computeDerivative gives the gradient of the wave function

    std::vector<double> qForce;

    for (int i=0; i < m_numberOfDimensions; i++){
        qForce.push_back(2*m_waveFunction->computeDerivative(m_particles)[particle*m_numberOfDimensions + i]);
    }
    return qForce;
}

double System::calculateGreensFunction(int particle, std::vector<double> positionOld, std::vector<double> positionNew){
    // Calculate the Greens function using the drift force from the quantumForce function

    double GreensFunction = 0;
    double D = 0.5;

    for (int i=0; i < m_numberOfDimensions; i++){
        double tmp = positionNew[i] - positionOld[i] - D*m_dt*quantumForce(particle)[i];
        GreensFunction += tmp*tmp;
    }

    GreensFunction = exp(-GreensFunction/4*D*m_dt);
    return GreensFunction;
}

void System::runMetropolisSteps(int numberOfMetropolisSteps, bool importanceSampling,
                                bool saveEnergies, bool savePositions, bool showProgress,
                                bool printToTerminal) {
    // Initialize Monte Carlo simulation
    m_particles                 = m_initialState->getParticles();
    m_sampler                   = new Sampler(this);
    m_numberOfMetropolisSteps   = numberOfMetropolisSteps;
    m_sampler->setNumberOfMetropolisSteps(numberOfMetropolisSteps);

    // Find CPU time
    clock_t start, finish;
    start = clock();
    int percent = numberOfMetropolisSteps/100;
    int progress = 0;
    double equilibrationSteps = m_equilibrationFraction*numberOfMetropolisSteps;

    for (int i=0; i < numberOfMetropolisSteps; i++) {

        // Update progress
        if (showProgress) {
            if (i%percent==0){
                progress += 1;
                cout << progress << "%" << "\n\033[F";
            }
        }

        bool acceptedStep;
        // Ask Metropolis step functions whether the step was accepted or not
        if (importanceSampling)  acceptedStep = metropolisStepImpSampling();
        else                     acceptedStep = metropolisStep();

        /* Here you should sample the energy (and maybe other things using
         * the m_sampler instance of the Sampler class. Make sure, though,
         * to only begin sampling after you have let the system equilibrate
         * for a while. You may handle this using the fraction of steps which
         * are equilibration steps; m_equilibrationFraction.
         */

        if (i>equilibrationSteps){
            // Sample energy etc.
            m_sampler->sample(acceptedStep, saveEnergies, savePositions);
        }
    }

    finish = clock();
    m_computationTime = (finish-start)/ (double) CLOCKS_PER_SEC;

    // Compute final expectation values
    m_sampler->computeAverages();
    if (printToTerminal) m_sampler->printOutputToTerminal();
}

void System::setNumberOfParticles(int numberOfParticles) {
    m_numberOfParticles = numberOfParticles;
}

void System::setNumberOfDimensions(int numberOfDimensions) {
    m_numberOfDimensions = numberOfDimensions;
}

void System::setStepLength(double stepLength) {
    assert(stepLength >= 0);
    m_stepLength = stepLength;
}

void System::setTimeStep(double dt) {
    assert(dt > 0);
    m_dt = dt;
}

void System::setEquilibrationFraction(double equilibrationFraction) {
    assert(equilibrationFraction >= 0);
    m_equilibrationFraction = equilibrationFraction;
}

void System::setHamiltonian(Hamiltonian* hamiltonian) {
    m_hamiltonian = hamiltonian;
}

void System::setWaveFunction(WaveFunction* waveFunction) {
    m_waveFunction = waveFunction;
}

void System::setInitialState(InitialState* initialState) {
    m_initialState = initialState;
}


