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

void SteepestDescent::obtainOptimalParameter(std::vector<double> parameters, double tol, int maxIterations,
                                             int numberOfMetropolisSteps, bool importanceSampling){

    int iteration = 0;  // Count iterations so we can force quit after a given number max iterations
    int numberOfParameters = parameters.size();
    double diff = 0;
    std::vector<double> parametersNew = parameters;
    //int parameterNumber;
    //if (parameterName == "alpha") parameterNumber = 0;
    //if (parameterName == "beta")  parameterNumber = 1;

    do{

        // Run Monte Carlo simulation to find expectation values
        m_system->getInitialState()->setupInitialState();
        for (int i=0; i < numberOfParameters; i++) {
            m_system->getWaveFunction()->adjustParameter(parameters[i], i);
        }

        m_system->runMetropolisSteps(numberOfMetropolisSteps, importanceSampling, false, false);

        std::vector<double> derivative(numberOfParameters);  //derivative of local energy.
        // Expectation values needed to calculate derivative of local energy:
        double energy = m_system->getSampler()->getEnergy();
        std::vector<double> waveFuncEnergy(numberOfParameters);
        std::vector<double> waveFuncDerivative(numberOfParameters);

        for (int i=0; i < numberOfParameters; i++) {
            waveFuncEnergy[i] = m_system->getSampler()->getWaveFuncEnergyParameters()[i];
            waveFuncDerivative[i] = m_system->getSampler()->getWaveFuncDerivativeParameters()[i];
            derivative[i] = 2*(waveFuncEnergy[i] - energy*waveFuncDerivative[i]);
        }

//        if (parameterNumber == 0){
//            waveFuncEnergy = m_system->getSampler()->getWaveFuncEnergyAlpha();
//            waveFuncDerivative = m_system->getSampler()->getWaveFuncDerivativeAlpha();
//        }

//        if (parameterNumber == 1){
//            waveFuncEnergy = m_system->getSampler()->getWaveFuncEnergyBeta();
//            waveFuncDerivative = m_system->getSampler()->getWaveFuncDerivativeBeta();
//        }

        // Find new parameter
        diff = 0;
        for (int i=0; i < numberOfParameters; i++) {
            parametersNew[i] = parameters[i] - derivative[i]*m_stepLengthSD;
            diff += abs(parametersNew[i] - parameters[i]);
        }
        //parametersNew = parameters - derivative*m_stepLengthSD;
        parameters = parametersNew;   // Update parameters
        iteration++;
        std::string upLine = "\e[A";
        cout << "Iterations: " << iteration << endl;
        for (int i=0; i < numberOfParameters; i++) {
            cout << "Parameter " << i+1 << ": " << parameters[i] << endl;
            upLine += "\e[A";
        }
        cout << upLine;             //"\033[F";



    }while(diff > tol && iteration < maxIterations);
    // Loop ends when requested tolerance for optimal parameters has been reached or after max iterations.

    cout << "Total iterations: " << iteration << endl;
    if (iteration==maxIterations) {
        cout << "Max iterations reached.\n";
        for (int i=0; i < numberOfParameters; i++) {
            cout << "Parameter" << i+1 << " at max iterations: " << parameters[i] << endl;
        }

    }
    else {
        for (int i=0; i < numberOfParameters; i++) {
            cout << "Optimal " << "Parameter" << i+1 << ": " << parameters[i] << endl;
        }
    }
    // Performing large MC simulation with optimal parameter:
    //m_system->getInitialState()->setupInitialState();
    for (int i=0; i < numberOfParameters; i++) {
        m_system->getWaveFunction()->adjustParameter(parameters[i], i);
    }
    //m_system->runMetropolisSteps((int) 1e6, importanceSampling, true, true);
}
