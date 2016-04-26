#include "manyelectrons.h"
#include <cmath>
#include <cassert>
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"
#include <iostream>

ManyElectrons::ManyElectrons(System* system, double alpha, double beta, double omega, double a, double C) :
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
    m_numberOfParticles = m_system->getNumberOfParticles();
    m_halfNumberOfParticles = m_numberOfParticles/2;
    setUpSlaterDet();
}

double ManyElectrons::evaluate(std::vector<class Particle*> particles) {


    //return waveFunction;
    return 0;
}

double ManyElectrons::evaluateSingleParticleWF(int nx, int ny, double x, double y) {

    double omegaSqrt = sqrt(m_omega);
    double waveFunction = computeHermitePolynomial(nx, x)*omegaSqrt*x
                         *computeHermitePolynomial(ny, y)*omegaSqrt*y
                         *exp(-m_omega*(x*x + y*y)*0.5);
    return waveFunction;
}

std::vector<double> ManyElectrons::computeDerivative(std::vector<class Particle*> particles){
    //Calculates ∇ψ/ψ for the interacting wave function using the analytical expression.


    //return derivative;
    //return 0;
}

double ManyElectrons::computeDoubleDerivative(std::vector<class Particle*> particles) {
    /* All wave functions need to implement this function, so you need to
     * find the double derivative analytically. Note that by double derivative,
     * we actually mean the sum of the Laplacians with respect to the
     * coordinates of each particle.
     *
     * This quantity is needed to compute the (local) energy (consider the
     * Schrödinger equation to see how the two are related).
     */

    //Calculates ∇²ψ/ψ for the interacting wave function using the analytical expression.

    //return doubleDerivative;
    return 0;
}

double ManyElectrons::computeDerivativeWrtAlpha(std::vector<Particle *> particles){
    // Calculates the derivative w.r.t. alpha for the interacting wave function using the analytical expression.


    return 0;
}

double ManyElectrons::computeDerivativeWrtBeta(std::vector<Particle *> particles){
    // Calculates the derivative w.r.t. beta for the interacting wave function using the analytical expression.


    return 0;
}

double ManyElectrons::computeMetropolisRatio(std::vector<Particle *> particles,
                                            int randomParticle, std::vector<double> positionChange) {

    int numberOfDimensions = m_system->getNumberOfDimensions();

    std::vector<double> positionOld = particles[randomParticle]->getPosition();

    for (int i=0; i<numberOfDimensions; i++) {
        particles[randomParticle]->adjustPosition(positionChange[i], i);
    }

    std::vector<double> positionNew = particles[randomParticle]->getPosition();

    int i = randomParticle;
    double ratioSlaterDet = 0;

    for (int j=0; j < m_halfNumberOfParticles; j++) {
        int nx = m_quantumNumbers(j, 0);
        int ny = m_quantumNumbers(j, 1);
        if (i < m_halfNumberOfParticles) {
            ratioSlaterDet += m_spinUpSlaterInverse(j,i)
                             *evaluateSingleParticleWF(nx, ny, positionNew[0], positionNew[1]);
        }
        else {
            ratioSlaterDet += m_spinDownSlaterInverse(j, i-m_halfNumberOfParticles)
                             *evaluateSingleParticleWF(nx, ny, positionNew[0], positionNew[1]);
        }
    }

    double beta = m_parameters[1];
    double exponent = 0;

    for (int j=0; j < m_numberOfParticles; j++) {
        if (j!=i) {
            double r_ijNew = 0;
            double r_ijOld = 0;
            for (int d=0; d < numberOfDimensions; d++) {
                r_ijNew += (positionNew[d] - particles[j]->getPosition()[d])
                          *(positionNew[d] - particles[j]->getPosition()[d]);
                r_ijOld += (positionOld[d] - particles[j]->getPosition()[d])
                          *(positionOld[d] - particles[j]->getPosition()[d]);
            }
            r_ijNew = sqrt(r_ijNew);
            r_ijOld = sqrt(r_ijOld);
            exponent += r_ijNew / (1 + beta*r_ijNew);
            exponent -= r_ijOld / (1 + beta*r_ijOld);
        }
    }
    double ratioJastrowFactor = exp(m_a*exponent);
    m_metropolisRatio = ratioSlaterDet*ratioJastrowFactor;
    return m_metropolisRatio;

}

double ManyElectrons::computeHermitePolynomial(int nValue, double position) {

    double HermitePolynomialPP = 0;                 // H_{n-2}
    double HermitePolynomialP = 1;                  // H_{n-1}
    double HermitePolynomial = HermitePolynomialP;  // H_n

    for (int n=1; n <= nValue; n++) {
        HermitePolynomial = 2*position*HermitePolynomialP - 2*(n-1)*HermitePolynomialPP;
        HermitePolynomialPP = HermitePolynomialP;
        HermitePolynomialP = HermitePolynomial;
    }

    return HermitePolynomial;

}

void ManyElectrons::setUpSlaterDet() {

    int n = 0;
    int nx = 0;
    int ny = 0;
    m_quantumNumbers = zeros<mat>(m_halfNumberOfParticles, 2);
    for (int p=0; p < m_halfNumberOfParticles; p++) {
        m_quantumNumbers(p, 0) = nx;    m_quantumNumbers(p, 1) = ny;
        if (ny == n) {
            n++;
            nx = n;
            ny = 0;
        }
        else {
            nx--;
            ny++;
        }
    }

    m_spinUpSlater = zeros<mat>(m_halfNumberOfParticles, m_halfNumberOfParticles);
    m_spinDownSlater = zeros<mat>(m_halfNumberOfParticles, m_halfNumberOfParticles);

    for (int i=0; i < m_halfNumberOfParticles; i++) {
        for (int j=0; j < m_halfNumberOfParticles; j++) {
            nx = m_quantumNumbers(j, 0);
            ny = m_quantumNumbers(j, 1);
            double xSpinUp = m_system->getParticles()[i]->getPosition()[0];
            double ySpinUp = m_system->getParticles()[i]->getPosition()[1];
            double xSpinDown = m_system->getParticles()[i+m_halfNumberOfParticles]->getPosition()[0];
            double ySpinDown = m_system->getParticles()[i+m_halfNumberOfParticles]->getPosition()[1];
            m_spinUpSlater(i,j) = evaluateSingleParticleWF(nx, ny, xSpinUp, ySpinUp);
            m_spinDownSlater(i,j) = evaluateSingleParticleWF(nx, ny, xSpinDown, ySpinDown);
        }
    }

    m_spinUpSlaterInverse = m_spinUpSlater.i();
    m_spinDownSlaterInverse = m_spinDownSlater.i();

}

void ManyElectrons::updateSlaterDet(int randomParticle) {
    int i = randomParticle;

    if (i < m_halfNumberOfParticles) {
        mat spinUpSlaterInverseOld = m_spinUpSlaterInverse;
        for (int j=0; j < m_halfNumberOfParticles; j++) {
            if (j!=i) {
                double sum = 0;
                for (int l=0; l <m_halfNumberOfParticles; l++) {
                    std::vector<double> r = m_system->getParticles()[l]->getPosition();
                    int nx = m_quantumNumbers(i, 0);
                    int ny = m_quantumNumbers(i, 1);
                    sum += evaluateSingleParticleWF(nx, ny, r[0], r[1])
                          *m_spinUpSlaterInverse(l,j);
                }
                for (int k=0; k < m_halfNumberOfParticles; k++) {
                    m_spinUpSlaterInverse(k,j) = spinUpSlaterInverseOld(k,j)
                                                -(sum/m_metropolisRatio)*spinUpSlaterInverseOld(k,i);
                }
            }
            else {
                for (int k=0; k < m_halfNumberOfParticles; k++) {
                    m_spinUpSlaterInverse(k,j) = spinUpSlaterInverseOld(k,i)/m_metropolisRatio;
                }
            }
        }
    }
    else {
        mat spinDownSlaterInverseOld = m_spinDownSlaterInverse;
        for (int j=0; j < m_halfNumberOfParticles; j++) {
            if (j != i-m_halfNumberOfParticles) {
                double sum = 0;
                for (int l=0; l < m_halfNumberOfParticles; l++) {
                    std::vector<double> r = m_system->getParticles()[l+m_halfNumberOfParticles]->getPosition();
                    int nx = m_quantumNumbers(i-m_halfNumberOfParticles, 0);
                    int ny = m_quantumNumbers(i-m_halfNumberOfParticles, 1);
                    sum += evaluateSingleParticleWF(nx, ny, r[0], r[1])
                          *m_spinDownSlaterInverse(l,j);
                }
                for (int k=0; k < m_halfNumberOfParticles; k++) {
                    m_spinDownSlaterInverse(k,j) = spinDownSlaterInverseOld(k,j)
                                                  -(sum/m_metropolisRatio)
                                                   *spinDownSlaterInverseOld(k, i-m_halfNumberOfParticles);
                }
            }
            else {
                for (int k=0; k < m_halfNumberOfParticles; k++) {
                    m_spinDownSlaterInverse(k,j) = spinDownSlaterInverseOld(k, i-m_halfNumberOfParticles)
                                                  /m_metropolisRatio;
                }
            }
        }
    }
}
