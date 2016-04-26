#include "twoelectrons.h"
#include <cmath>
#include <cassert>
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"
#include <iostream>

TwoElectrons::TwoElectrons(System* system, double alpha, double beta, double omega, double a, double C) :
        WaveFunction(system) {
    assert(omega > 0);
    m_omega = omega;
    assert(alpha >= 0);
    m_numberOfParameters = 2;
    m_parameters.reserve(2);
    m_parameters.push_back(alpha);
    m_parameters.push_back(beta);
    m_a = a;
    m_C = C;
}

double TwoElectrons::evaluate(std::vector<class Particle*> particles) {

    int numberOfParticles = m_system->getNumberOfParticles();
    int numberOfDimensions = m_system->getNumberOfDimensions();
    assert(numberOfParticles = 2);
    //assert(numberOfDimensions) = 2;
    double alpha = m_parameters[0];
    double beta = m_parameters[1];

    double r1Squared = 0;
    double r2Squared = 0;
    double r12 = 0;

    for (int i=0; i < numberOfDimensions; i++){
        std::vector<double> r1 = particles[0]->getPosition();
        std::vector<double> r2 = particles[1]->getPosition();

        r1Squared += r1[i]*r1[i];
        r2Squared += r2[i]*r2[i];
        r12 += (r1[i]-r2[i])*(r1[i]-r2[i]);
    }

    r12 = sqrt(r12);
    double waveFunction = m_C*exp(-alpha*m_omega*(r1Squared+r2Squared)*0.5)*exp(m_a*r12/(1+beta*r12));

    return waveFunction;
    //return 0;
}

std::vector<double> TwoElectrons::computeDerivative(std::vector<class Particle*> particles){
    //Calculates ∇ψ/ψ for the interacting wave function using the analytical expression.

    int numberOfParticles = m_system->getNumberOfParticles();
    int numberOfDimensions = m_system->getNumberOfDimensions();
    assert(numberOfParticles = 2);
    //assert(numberOfDimensions) = 2;
    double alpha = m_parameters[0];
    double beta = m_parameters[1];

    std::vector<double> r1 = particles[0]->getPosition();
    std::vector<double> r2 = particles[1]->getPosition();
    double r12 = 0;

    for (int i=0; i < numberOfDimensions; i++){
        r12 += (r1[i]-r2[i])*(r1[i]-r2[i]);
    }

    std::vector<double> derivative(numberOfParticles*numberOfDimensions);
    double constFactor = m_a/(r12*(1+beta*r12)*(1+beta*r12));

    for (int i=0; i < numberOfParticles; i++){
        for (int j=0; j < numberOfDimensions; j++){
            if (i==0) derivative[i+j] = -alpha*m_omega*r1[j] + (r1[j]-r2[j])*constFactor;
            else      derivative[i+j] = -alpha*m_omega*r1[j] - (r1[j]-r2[j])*constFactor;
        }
    }

    return derivative;
}

double TwoElectrons::computeDoubleDerivative(std::vector<class Particle*> particles) {
    /* All wave functions need to implement this function, so you need to
     * find the double derivative analytically. Note that by double derivative,
     * we actually mean the sum of the Laplacians with respect to the
     * coordinates of each particle.
     *
     * This quantity is needed to compute the (local) energy (consider the
     * Schrödinger equation to see how the two are related).
     */

    //Calculates ∇²ψ/ψ for the interacting wave function using the analytical expression.
    int numberOfParticles = m_system->getNumberOfParticles();
    int numberOfDimensions = m_system->getNumberOfDimensions();
    assert(numberOfParticles = 2);
    //assert(numberOfDimensions) = 2;
    double alpha = m_parameters[0];
    double beta = m_parameters[1];

    double r1Squared = 0;
    double r2Squared = 0;
    double r12 = 0;

    for (int i=0; i < numberOfDimensions; i++){
        std::vector<double> r1 = particles[0]->getPosition();
        std::vector<double> r2 = particles[1]->getPosition();

        r1Squared += r1[i]*r1[i];
        r2Squared += r2[i]*r2[i];
        r12 += (r1[i]-r2[i])*(r1[i]-r2[i]);
    }

    r12 = sqrt(r12);
    double denom = 1+beta*r12;

    double doubleDerivative = -2*numberOfDimensions*alpha*m_omega
                            + alpha*alpha*m_omega*m_omega*(r1Squared+r2Squared)
                            + 2*m_a*m_a/(denom*denom*denom*denom)
                            - alpha*m_omega*r12*m_a/(denom*denom);

    return doubleDerivative;
    //return 0;
}

double TwoElectrons::computeDerivativeWrtAlpha(std::vector<Particle *> particles){
    // Calculates the derivative w.r.t. alpha for the interacting wave function using the analytical expression.

    int numberOfParticles = m_system->getNumberOfParticles();
    int numberOfDimensions = m_system->getNumberOfDimensions();
    assert(numberOfParticles = 2);
    //assert(numberOfDimensions) = 2;
    double alpha = m_parameters[0];
    double beta = m_parameters[1];

    double r1Squared = 0;
    double r2Squared = 0;
    double r12 = 0;

    for (int i=0; i < numberOfDimensions; i++){
        std::vector<double> r1 = particles[0]->getPosition();
        std::vector<double> r2 = particles[1]->getPosition();

        r1Squared += r1[i]*r1[i];
        r2Squared += r2[i]*r2[i];
        r12 += (r1[i]-r2[i])*(r1[i]-r2[i]);
    }

    r12 = sqrt(r12);

    double derivative = -m_C*m_omega*(r1Squared + r2Squared)*0.5
                        *exp(-alpha*m_omega*(r1Squared + r2Squared)*0.5)
                        *exp(m_a*r12/(1+beta*r12));

    return derivative;
}

double TwoElectrons::computeDerivativeWrtBeta(std::vector<Particle *> particles){
    // Calculates the derivative w.r.t. beta for the interacting wave function using the analytical expression.

    int numberOfParticles = m_system->getNumberOfParticles();
    int numberOfDimensions = m_system->getNumberOfDimensions();
    assert(numberOfParticles = 2);
    //assert(numberOfDimensions) = 2;
    double alpha = m_parameters[0];
    double beta = m_parameters[1];

    double r1Squared = 0;
    double r2Squared = 0;
    double r12 = 0;

    for (int i=0; i < numberOfDimensions; i++){
        std::vector<double> r1 = particles[0]->getPosition();
        std::vector<double> r2 = particles[1]->getPosition();

        r1Squared += r1[i]*r1[i];
        r2Squared += r2[i]*r2[i];
        r12 += (r1[i]-r2[i])*(r1[i]-r2[i]);
    }

    r12 = sqrt(r12);

    double derivative = -m_C*m_a*r12*r12/((1+beta*r12)*(1+beta*r12))
                        *exp(-alpha*m_omega*(r1Squared + r2Squared)*0.5)
                        *exp(m_a*r12/(1+beta*r12));

    return derivative;
}

double TwoElectrons::computeMetropolisRatio(std::vector<Particle *> particles,
                                            int randomParticle, std::vector<double> positionChange) {
    int numberOfDimensions = m_system->getNumberOfDimensions();

    // Evaluate the wave function for current positions
    double waveFunctionOld = evaluate(particles);

    // Change position to trial state
    for (int i=0; i<numberOfDimensions; i++){
        particles[randomParticle]->adjustPosition(positionChange[i], i);
    }

    // Evaluate the wave function for the trial state
    double waveFunctionNew = evaluate(particles);

    return waveFunctionNew*waveFunctionNew / (waveFunctionOld*waveFunctionOld);
}
