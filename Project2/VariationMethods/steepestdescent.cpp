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

void SteepestDescent::obtainOptimalParameter(double parameter, std::string parameterName,
                                         double tol, int maxIterations, int numberOfMetropolisSteps,
                                         bool importanceSampling){

    int iteration = 0;  // Count iterations so we can force quit after a given number max iterations
    double parameterNew = parameter;
    int parameterNumber;
    if (parameterName == "alpha") parameterNumber = 0;
    if (parameterName == "beta")  parameterNumber = 1;

    do{
        parameter = parameterNew;   // Update alpha

        // Run Monte Carlo simulation to find expectation values
        m_system->getInitialState()->setupInitialState();
        m_system->getWaveFunction()->adjustParameter(parameter, parameterNumber);
        m_system->runMetropolisSteps(numberOfMetropolisSteps, importanceSampling, false, false);

        double derivative;  //derivative of local energy.
        // Expectation values needed to calculate derivative of local energy:
        double energy = m_system->getSampler()->getEnergy();
        double waveFuncEnergy;
        double waveFuncDerivative;

        if (parameterNumber == 0){
            waveFuncEnergy = m_system->getSampler()->getWaveFuncEnergyAlpha();
            waveFuncDerivative = m_system->getSampler()->getWaveFuncDerivativeAlpha();
        }

        if (parameterNumber == 1){
            waveFuncEnergy = m_system->getSampler()->getWaveFuncEnergyBeta();
            waveFuncDerivative = m_system->getSampler()->getWaveFuncDerivativeBeta();
        }

        derivative = 2*(waveFuncEnergy - energy*waveFuncDerivative);
        // Find new parameter
        parameterNew = parameter - derivative*m_stepLengthSD;
        iteration++;
        cout << "Iterations: " << iteration << endl;
        cout << parameterName << ": " << parameter << "\033[F";

    }while(abs(parameterNew - parameter) > tol && iteration < maxIterations);
    // Loop ends when requested tolerance for optimal alpha has been reached or after max iterations.

    cout << "Total iterations: " << iteration << endl;
    //const char* message = "Optimal " << parameterName << ": ";
    if (iteration==maxIterations) cout << "Max iterations reached.\n"<< parameterName
                                       << " at max iterations: " << parameter << endl;
    else cout << "Optimal " << parameterName << ": " << parameter << endl;

    // Performing large MC simulation with optimal parameter:
    //m_system->getInitialState()->setupInitialState();
    m_system->getWaveFunction()->adjustParameter(parameter, parameterNumber);
    //m_system->runMetropolisSteps((int) 1e6, importanceSampling, true, true);
}
