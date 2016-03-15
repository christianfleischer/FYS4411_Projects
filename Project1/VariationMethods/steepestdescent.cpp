#include "steepestdescent.h"
#include "system.h"
#include "sampler.h"
#include "InitialStates/randomuniform.h"
#include "WaveFunctions/wavefunction.h"
#include <cmath>
#include <iostream>

using namespace std;

SteepestDescent::SteepestDescent(System* system, double stepLengthSD)
{
    m_system = system;
    m_stepLengthSD = stepLengthSD;
}

void SteepestDescent::obtainOptimalAlpha(double alpha, double tol,
                                         int maxIterations, int numberOfMetropolisSteps,
                                         bool importanceSampling){

    int iteration = 0;  // Count iterations so we can force quit after a given number max iterations
    double alphaNew = alpha;

    do{
        alpha = alphaNew;   // Update alpha

        // Run Monte Carlo simulation to find expectation values
        m_system->getInitialState()->setupInitialState();
        m_system->getWaveFunction()->adjustParameter(alpha, 0);
        m_system->runMetropolisSteps(numberOfMetropolisSteps, importanceSampling, false, false, false, false);

        double derivative;  //derivative of local energy.
        // Expectation values needed to calculate derivative of local energy:
        double energy = m_system->getSampler()->getEnergy();
        double waveFuncEnergy = m_system->getSampler()->getWaveFuncEnergy();
        double waveFuncDerivative = m_system->getSampler()->getWaveFuncDerivative();

        derivative = 2*(waveFuncEnergy - energy*waveFuncDerivative);
        // Find new alpha
        alphaNew = alpha - derivative*m_stepLengthSD;
        iteration++;
        cout << "Iterations: " << iteration << endl;
        cout << "Alpha: " << alpha << "\033[F";

    }while(abs(alphaNew - alpha) > tol && iteration < maxIterations);
    // Loop ends when requested tolerance for optimal alpha has been reached or after max iterations.

    cout << "Total iterations: " << iteration << endl;
    const char* message = "Optimal alpha: ";
    if (iteration==maxIterations) message = "Max iterations reached.\nAlpha at max iterations: ";
    cout << message << alpha << endl;

    // Performing large MC simulation with optimal aplha:
    m_system->getInitialState()->setupInitialState();
    m_system->getWaveFunction()->adjustParameter(alpha, 0);
    m_system->runMetropolisSteps((int) 1e6, importanceSampling, false, false, true, true);
}
