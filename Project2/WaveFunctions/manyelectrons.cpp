#include "manyelectrons.h"
#include <cmath>
#include <cassert>
#include "../InitialStates/randomuniform.h"
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"
#include <iostream>

ManyElectrons::ManyElectrons(System* system, double alpha, double beta, double omega, double C, bool Jastrow) :
        WaveFunction(system) {
    assert(omega > 0);
    m_omega = omega;
    assert(alpha >= 0);
    m_Jastrow = Jastrow;
    m_numberOfParameters = 2;
    m_parameters.reserve(2);
    m_parameters.push_back(alpha);
    m_parameters.push_back(beta);
    m_C = C;
    m_numberOfParticles = m_system->getNumberOfParticles();
    m_halfNumberOfParticles = m_numberOfParticles/2.;
    setUpSlaterDet();
}

double ManyElectrons::evaluate(std::vector<class Particle*> particles) {
    // Evaluates the wave function using brute force.
    mat spinUpSlater = zeros<mat>(m_halfNumberOfParticles, m_halfNumberOfParticles);
    mat spinDownSlater = zeros<mat>(m_halfNumberOfParticles, m_halfNumberOfParticles);
    double beta = m_parameters[1];

    for (int i=0; i < m_halfNumberOfParticles; i++) {
        for (int j=0; j < m_halfNumberOfParticles; j++) {
            int nx = m_quantumNumbers(j, 0);
            int ny = m_quantumNumbers(j, 1);
            double xSpinUp = m_system->getParticles()[i]->getPosition()[0];
            double ySpinUp = m_system->getParticles()[i]->getPosition()[1];
            double xSpinDown = m_system->getParticles()[i+m_halfNumberOfParticles]->getPosition()[0];
            double ySpinDown = m_system->getParticles()[i+m_halfNumberOfParticles]->getPosition()[1];
            spinUpSlater(i,j) = evaluateSingleParticleWF(nx, ny, xSpinUp, ySpinUp);
            spinDownSlater(i,j) = evaluateSingleParticleWF(nx, ny, xSpinDown, ySpinDown);
        }
    }

    double exponent = 0;
    if (m_Jastrow) {
        for (int i=0; i < m_numberOfParticles; i++) {
            std::vector<double> r_i = m_system->getParticles()[i]->getPosition();
            for (int j=i+1; j < m_numberOfParticles; j++) {
                std::vector<double> r_j = m_system->getParticles()[j]->getPosition();
                double r_ij = (r_i[0] - r_j[0])*(r_i[0] - r_j[0]) + (r_i[1] - r_j[1])*(r_i[1] - r_j[1]);
                r_ij = sqrt(r_ij);
                double denom = 1+beta*r_ij;
                exponent += m_a(i,j)*r_ij/denom;
            }
        }
    }

    double waveFunction = det(spinDownSlater)*det(spinUpSlater)*exp(exponent);

    return waveFunction;
}

double ManyElectrons::evaluateSingleParticleWF(int nx, int ny, double x, double y) {
    // Calculates the single particle wave function.
    double alpha = m_parameters[0];

    double waveFunction = computeHermitePolynomial(nx, x)
                         *computeHermitePolynomial(ny, y)
                         *exp(-m_omega*alpha*(x*x + y*y)*0.5);
    return waveFunction;
}

std::vector<double> ManyElectrons::computeDerivative(std::vector<class Particle*> particles,
                                                     int randomParticle) {
    //Calculates ∇ψ/ψ for the wave function.

    int i = randomParticle;
    int numberOfParticles = m_system->getNumberOfParticles();
    int numberOfDimensions = m_system->getNumberOfDimensions();
    std::vector<double> derivative(numberOfParticles*numberOfDimensions);
    derivative[i*numberOfDimensions] = computeSlaterGradient(particles, i)[0]
                                        ;//+computeJastrowGradient(particles, i)[0];
    derivative[i*numberOfDimensions+1] = computeSlaterGradient(particles, i)[1]
                                          ;//+computeJastrowGradient(particles, i)[1];
    if (m_Jastrow) {
        derivative[i*numberOfDimensions] += computeJastrowGradient(particles, i)[0];
        derivative[i*numberOfDimensions+1] += computeJastrowGradient(particles, i)[1];
    }
    return derivative;
    //return 0;
}

std::vector<double> ManyElectrons::computeSlaterGradient(std::vector<class Particle*> particles, int i) {
    // Computes the gradient of the Slater part of the wave function.
    std::vector<double> slaterGradient(2);
    slaterGradient[0] = 0;
    slaterGradient[1] = 0;
    double x = particles[i]->getPosition()[0];
    double y = particles[i]->getPosition()[1];

    if (i < m_halfNumberOfParticles) {
        for (int j=0; j < m_halfNumberOfParticles; j++) {
            int nx = m_quantumNumbers(j, 0);
            int ny = m_quantumNumbers(j, 1);
            std::vector<double> SPWFGradient = computeSPWFDerivative(nx, ny, x, y);
            slaterGradient[0] += SPWFGradient[0]*m_spinUpSlaterInverse(j,i);
            slaterGradient[1] += SPWFGradient[1]*m_spinUpSlaterInverse(j,i);
        }
    }
    else {
        for (int j=0; j < m_halfNumberOfParticles; j++) {
            int nx = m_quantumNumbers(j, 0);
            int ny = m_quantumNumbers(j, 1);
            std::vector<double> SPWFGradient = computeSPWFDerivative(nx, ny, x, y);
            slaterGradient[0] += SPWFGradient[0]*m_spinDownSlaterInverse(j, i-m_halfNumberOfParticles);
            slaterGradient[1] += SPWFGradient[1]*m_spinDownSlaterInverse(j, i-m_halfNumberOfParticles);
        }
    }

    return slaterGradient;

}

std::vector<double> ManyElectrons::computeJastrowGradient(std::vector<class Particle*> particles, int k) {
    // Computes the gradient of the Jastrow part of the wave function.
    std::vector<double> jastrowGradient(2);
    jastrowGradient[0] = jastrowGradient[1] = 0;

    double beta = m_parameters[1];

    for (int j=0; j < k; j++) {
        std::vector<double> r_k = particles[k]->getPosition();
        std::vector<double> r_j = particles[j]->getPosition();
        double r_kj = (r_k[0]-r_j[0])*(r_k[0]-r_j[0]) + (r_k[1]-r_j[1])*(r_k[1]-r_j[1]);
        r_kj = sqrt(r_kj);
        double denom = 1 + beta*r_kj;
        jastrowGradient[0] += (r_k[0]-r_j[0])/r_kj * m_a(k, j)/(denom*denom);
        jastrowGradient[1] += (r_k[1]-r_j[1])/r_kj * m_a(k, j)/(denom*denom);
    }

    for (int j=k+1; j < m_numberOfParticles; j++) {
        std::vector<double> r_k = particles[k]->getPosition();
        std::vector<double> r_j = particles[j]->getPosition();
        double r_kj = (r_k[0]-r_j[0])*(r_k[0]-r_j[0]) + (r_k[1]-r_j[1])*(r_k[1]-r_j[1]);
        r_kj = sqrt(r_kj);
        double denom = 1 + beta*r_kj;
        jastrowGradient[0] += (r_k[0]-r_j[0])/r_kj * m_a(k, j)/(denom*denom);
        jastrowGradient[1] += (r_k[1]-r_j[1])/r_kj * m_a(k, j)/(denom*denom);
    }

    return jastrowGradient;

}

std::vector<double> ManyElectrons::computeSPWFDerivative(int nx, int ny, double x, double y) {
    // Calculates the single particle wave function differentiated w.r.t. position.
    std::vector<double> derivative(m_system->getNumberOfDimensions());
    double alpha = m_parameters[0];
    double r2 = x*x + y*y;

    derivative[0] = (computeHermitePolynomialDerivative(nx, x) - alpha*m_omega*x*computeHermitePolynomial(nx, x))
                   *computeHermitePolynomial(ny, y)*exp(-alpha*m_omega*r2*0.5);

    derivative[1] = (computeHermitePolynomialDerivative(ny, y) - alpha*m_omega*y*computeHermitePolynomial(ny, y))
                   *computeHermitePolynomial(nx, x)*exp(-alpha*m_omega*r2*0.5);

    return derivative;
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

    //Calculates ∇²ψ/ψ for the wave function.

    int numberOfDimensions = m_system->getNumberOfDimensions();
    double slaterLaplacian = 0;
    double jastrowLaplacian = 0;
    double crossTerm = 0;

    for (int i=0; i < m_halfNumberOfParticles; i++) {
        for (int j=0; j < m_halfNumberOfParticles; j++) {
            int nx = m_quantumNumbers(j, 0);
            int ny = m_quantumNumbers(j, 1);
            double x = particles[i]->getPosition()[0];
            double y = particles[i]->getPosition()[1];
            slaterLaplacian += computeSPWFDoubleDerivative(nx, ny, x, y)*m_spinUpSlaterInverse(j,i);
        }
    }
    for (int i=m_halfNumberOfParticles; i < m_numberOfParticles; i++) {
        for (int j=0; j < m_halfNumberOfParticles; j++) {
            int nx = m_quantumNumbers(j, 0);
            int ny = m_quantumNumbers(j, 1);
            double x = particles[i]->getPosition()[0];
            double y = particles[i]->getPosition()[1];
            slaterLaplacian += computeSPWFDoubleDerivative(nx, ny, x, y)
                              *m_spinDownSlaterInverse(j,i-m_halfNumberOfParticles);
        }
    }

    if (m_Jastrow) {
        double beta = m_parameters[1];
        double jastrowSum1 = 0;
        double jastrowSum2 = 0;
        int d = numberOfDimensions;

        for (int k=0; k < m_numberOfParticles; k++) {
            std::vector<double> r_k = particles[k]->getPosition();
            for (int j=0; j < k; j++) {
                std::vector<double> r_j = particles[j]->getPosition();
                double r_kj = (r_k[0]-r_j[0])*(r_k[0]-r_j[0]) + (r_k[1]-r_j[1])*(r_k[1]-r_j[1]);
                r_kj = sqrt(r_kj);
                double denom_kj = (1+beta*r_kj);
                jastrowSum2 += (d-1)*m_a(k,j)/(r_kj*denom_kj*denom_kj) - 2*m_a(k,j)*beta/(denom_kj*denom_kj*denom_kj);

                for (int i=0; i < k; i++) {
                    std::vector<double> r_i = particles[i]->getPosition();
                    double r_ki = (r_k[0]-r_i[0])*(r_k[0]-r_i[0]) + (r_k[1]-r_i[1])*(r_k[1]-r_i[1]);
                    r_ki = sqrt(r_ki);
                    double denom_ki = (1+beta*r_ki);
                    double factor1 = (r_k[0]-r_i[0])*(r_k[0]-r_j[0]) + (r_k[1]-r_i[1])*(r_k[1]-r_j[1]);
                    factor1 /= r_ki*r_kj;
                    double factor2 = m_a(k,i)/(denom_ki*denom_ki) * m_a(k,j)/(denom_kj*denom_kj);
                    jastrowSum1 += factor1*factor2;
                }
                for (int i=k+1; i < m_numberOfParticles; i++) {
                    std::vector<double> r_i = particles[i]->getPosition();
                    double r_ki = (r_k[0]-r_i[0])*(r_k[0]-r_i[0]) + (r_k[1]-r_i[1])*(r_k[1]-r_i[1]);
                    r_ki = sqrt(r_ki);
                    double denom_ki = (1+beta*r_ki);
                    double factor1 = (r_k[0]-r_i[0])*(r_k[0]-r_j[0]) + (r_k[1]-r_i[1])*(r_k[1]-r_j[1]);
                    factor1 /= r_ki*r_kj;
                    double factor2 = m_a(k,i)/(denom_ki*denom_ki) * m_a(k,j)/(denom_kj*denom_kj);
                    jastrowSum1 += factor1*factor2;
                }
            }
            for (int j=k+1; j < m_numberOfParticles; j++) {
                std::vector<double> r_j = particles[j]->getPosition();
                double r_kj = (r_k[0]-r_j[0])*(r_k[0]-r_j[0]) + (r_k[1]-r_j[1])*(r_k[1]-r_j[1]);
                r_kj = sqrt(r_kj);
                double denom_kj = (1+beta*r_kj);
                jastrowSum2 += (d-1)*m_a(k,j)/(r_kj*denom_kj*denom_kj) - 2*m_a(k,j)*beta/(denom_kj*denom_kj*denom_kj);

                for (int i=0; i < k; i++) {
                    std::vector<double> r_i = particles[i]->getPosition();
                    double r_ki = (r_k[0]-r_i[0])*(r_k[0]-r_i[0]) + (r_k[1]-r_i[1])*(r_k[1]-r_i[1]);
                    r_ki = sqrt(r_ki);
                    double denom_ki = (1+beta*r_ki);
                    double factor1 = (r_k[0]-r_i[0])*(r_k[0]-r_j[0]) + (r_k[1]-r_i[1])*(r_k[1]-r_j[1]);
                    factor1 /= r_ki*r_kj;
                    double factor2 = m_a(k,i)/(denom_ki*denom_ki) * m_a(k,j)/(denom_kj*denom_kj);
                    jastrowSum1 += factor1*factor2;
                }
                for (int i=k+1; i < m_numberOfParticles; i++) {
                    std::vector<double> r_i = particles[i]->getPosition();
                    double r_ki = (r_k[0]-r_i[0])*(r_k[0]-r_i[0]) + (r_k[1]-r_i[1])*(r_k[1]-r_i[1]);
                    r_ki = sqrt(r_ki);
                    double denom_ki = (1+beta*r_ki);
                    double factor1 = (r_k[0]-r_i[0])*(r_k[0]-r_j[0]) + (r_k[1]-r_i[1])*(r_k[1]-r_j[1]);
                    factor1 /= r_ki*r_kj;
                    double factor2 = m_a(k,i)/(denom_ki*denom_ki) * m_a(k,j)/(denom_kj*denom_kj);
                    jastrowSum1 += factor1*factor2;
                }
            }
        }

        jastrowLaplacian = jastrowSum1 + jastrowSum2;

        for (int d=0; d < numberOfDimensions; d++) {
            for (int i=0; i < m_numberOfParticles; i++) {
                crossTerm += computeSlaterGradient(particles, i)[d]*computeJastrowGradient(particles, i)[d];
            }
        }
    }

    double laplacian = slaterLaplacian + jastrowLaplacian + 2*crossTerm;
    return laplacian;
    //return 0;
}

double ManyElectrons::computeSPWFDoubleDerivative(int nx, int ny, double x, double y) {
    // Calculates the single particle wave function twice differentiated w.r.t. position.
    double doubleDerivative = 0;
    double alpha = m_parameters[0];
    double r2 = x*x + y*y;

    doubleDerivative += computeHermitePolynomial(ny, y)*exp(-alpha*m_omega*r2*0.5)
                       *(computeHermitePolynomialDoubleDerivative(nx, x)
                         - alpha*m_omega*computeHermitePolynomial(nx, x)
                         - 2*alpha*m_omega*x*computeHermitePolynomialDerivative(nx, x)
                         + alpha*alpha*m_omega*m_omega*x*x*computeHermitePolynomial(nx, x));

    doubleDerivative += computeHermitePolynomial(nx, x)*exp(-alpha*m_omega*r2*0.5)
                       *(computeHermitePolynomialDoubleDerivative(ny, y)
                         - alpha*m_omega*computeHermitePolynomial(ny, y)
                         - 2*alpha*m_omega*y*computeHermitePolynomialDerivative(ny, y)
                         + alpha*alpha*m_omega*m_omega*y*y*computeHermitePolynomial(ny, y));

    return doubleDerivative;

}

double ManyElectrons::computeSPWFAlphaDerivative(int nx, int ny, double x, double y) {
    // Calculates the single particle wave function differentiated w.r.t. alpha.
    double derivative = 0;
    double alpha = m_parameters[0];
    double r2 = x*x + y*y;

    derivative += (-0.5*m_omega*r2*computeHermitePolynomial(nx, x)*computeHermitePolynomial(ny, y)
                   +computeHermitePolynomialAlphaDerivative(nx, x)*computeHermitePolynomial(ny, y)
                   +computeHermitePolynomialAlphaDerivative(ny, y)*computeHermitePolynomial(nx, x))
                 *exp(-0.5*alpha*m_omega*r2);

    return derivative;

}

std::vector<double> ManyElectrons::computeDerivativeWrtParameters(std::vector<Particle *> particles){
    // Calculates the derivative w.r.t. alpha and beta for the interacting wave function using the analytical expression.

    std::vector<double> derivative(2);
    double slaterUpAlphaDerivative = 0;
    double slaterDownAlphaDerivative = 0;

    for (int i=0; i < m_halfNumberOfParticles; i++) {
        for (int j=0; j < m_halfNumberOfParticles; j++) {
            int nx = m_quantumNumbers(j, 0);
            int ny = m_quantumNumbers(j, 1);
            double xSpinUp = particles[i]->getPosition()[0];
            double ySpinUp = particles[i]->getPosition()[1];
            double xSpinDown = particles[i+m_halfNumberOfParticles]->getPosition()[0];
            double ySpinDown = particles[i+m_halfNumberOfParticles]->getPosition()[1];
            slaterUpAlphaDerivative += computeSPWFAlphaDerivative(nx, ny, xSpinUp, ySpinUp)
                                       *m_spinUpSlaterInverse(j,i);
            slaterDownAlphaDerivative += computeSPWFAlphaDerivative(nx, ny, xSpinDown, ySpinDown)
                                         *m_spinDownSlaterInverse(j,i);
        }
    }

    double beta = m_parameters[1];
    int numberOfDimensions = m_system->getNumberOfDimensions();
    double exponent = 0;
    double betaDerivative = 0;
    double r_ij = 0;

    for (int i=0; i < m_numberOfParticles; i++) {
        std::vector<double> r_i = particles[i]->getPosition();
        for (int j=i+1; j < m_numberOfParticles; j++) {
            std::vector<double> r_j = particles[j]->getPosition();
            for (int k=0; k < numberOfDimensions; k++) {
                r_ij += (r_i[k]-r_j[k])*(r_i[k]-r_j[k]);
            }
            r_ij = sqrt(r_ij);

            exponent += m_a(i,j)*r_ij/(1. + beta*r_ij);
            double denom = (1+beta*r_ij);
            betaDerivative -= m_a(i,j)*r_ij*r_ij/(denom*denom);
        }
    }

    derivative[0] = (slaterUpAlphaDerivative + slaterDownAlphaDerivative)*evaluate(particles);//*exp(exponent);
    derivative[1] = betaDerivative*evaluate(particles);

    return derivative;
    //return 0;
}

double ManyElectrons::computeMetropolisRatio(std::vector<Particle *> particles,
                                            int randomParticle, std::vector<double> positionChange) {
    // Function for calculating the wave function part of the Metropolis ratio,
    // both the Slater part and the Jastrow part.
    int numberOfDimensions = m_system->getNumberOfDimensions();

    std::vector<double> positionOld = particles[randomParticle]->getPosition();

    for (int i=0; i<numberOfDimensions; i++) {
        particles[randomParticle]->adjustPosition(positionChange[i], i);
    }

    std::vector<double> positionNew = particles[randomParticle]->getPosition();

    int i = randomParticle;
    double ratioSlaterDet = 0;

    if (i < m_halfNumberOfParticles) {
        for (int j=0; j < m_halfNumberOfParticles; j++) {
            int nx = m_quantumNumbers(j, 0);
            int ny = m_quantumNumbers(j, 1);
            ratioSlaterDet += m_spinUpSlaterInverse(j,i)
                             *evaluateSingleParticleWF(nx, ny, positionNew[0], positionNew[1]);
        }
    }
    else {
        for (int j=0; j < m_halfNumberOfParticles; j++) {
            int nx = m_quantumNumbers(j, 0);
            int ny = m_quantumNumbers(j, 1);
            ratioSlaterDet += m_spinDownSlaterInverse(j, i-m_halfNumberOfParticles)
                             *evaluateSingleParticleWF(nx, ny, positionNew[0], positionNew[1]);
        }
    }

    double beta = m_parameters[1];
    double exponent = 0;

    if (m_Jastrow) {
        for (int j=0; j < i; j++) {
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

            exponent += m_a(i,j)*r_ijNew / (1. + beta*r_ijNew);
            exponent -= m_a(i,j)*r_ijOld / (1. + beta*r_ijOld);
        }
        for (int j=i+1; j < m_numberOfParticles; j++) {
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

            exponent += m_a(i,j)*r_ijNew / (1. + beta*r_ijNew);
            exponent -= m_a(i,j)*r_ijOld / (1. + beta*r_ijOld);
        }
    }
    double ratioJastrowFactor = exp(exponent);
    m_ratioSlaterDet = ratioSlaterDet;
    m_metropolisRatio = ratioSlaterDet*ratioJastrowFactor;

    return m_metropolisRatio;

}

double ManyElectrons::computeHermitePolynomial(int nValue, double position) {
    // Computes Hermite polynomials.
    double alpha = m_parameters[0];
    double alphaSqrt = sqrt(alpha);
    double omegaSqrt = sqrt(m_omega);
    double variable = alphaSqrt*omegaSqrt*position;

    double HermitePolynomialPP = 0;                 // H_{n-2}
    double HermitePolynomialP = 1;                  // H_{n-1}
    double HermitePolynomial = HermitePolynomialP;  // H_n

    for (int n=1; n <= nValue; n++) {
        HermitePolynomial = 2*variable*HermitePolynomialP - 2*(n-1)*HermitePolynomialPP;
        HermitePolynomialPP = HermitePolynomialP;
        HermitePolynomialP = HermitePolynomial;
    }

    return HermitePolynomial;

}

double ManyElectrons::computeHermitePolynomialDerivative(int nValue, double position) {
    // Computes Hermite polynomials differentiated w.r.t. position.
    double alpha = m_parameters[0];
    double alphaSqrt = sqrt(alpha);
    double omegaSqrt = sqrt(m_omega);

    double HPDerivativePP = 0;              // d/dx H_{n-2}
    double HPDerivativeP = 0;               // d/dx H_{n-1}
    double HPDerivative = HPDerivativeP;    // d/dx H_n

    for (int n=1; n <= nValue; n++) {
        HPDerivative = 2*alphaSqrt*omegaSqrt*computeHermitePolynomial(n-1, position)
                      +2*alphaSqrt*omegaSqrt*position*HPDerivativeP
                      -2*(n-1)*HPDerivativePP;
        HPDerivativePP = HPDerivativeP;
        HPDerivativeP = HPDerivative;
    }

    return HPDerivative;

}

double ManyElectrons::computeHermitePolynomialDoubleDerivative(int nValue, double position) {
    // Computes Hermite polynomials twice differentiated w.r.t. position.
    double alpha = m_parameters[0];
    double alphaSqrt = sqrt(alpha);
    double omegaSqrt = sqrt(m_omega);

    double HPDoubleDerivativePP = 0;                    // d/dx d/dx H_{n-2}
    double HPDoubleDerivativeP = 0;                     // d/dx d/dx H_{n-1}
    double HPDoubleDerivative = HPDoubleDerivativeP;    // d/dx d/dx H_n

    for (int n=1; n <= nValue; n++) {
        HPDoubleDerivative = 4*alphaSqrt*omegaSqrt*computeHermitePolynomialDerivative(n-1, position)
                            +2*alphaSqrt*omegaSqrt*position*HPDoubleDerivativeP
                            -2*(n-1)*HPDoubleDerivativePP;
        HPDoubleDerivativePP = HPDoubleDerivativeP;
        HPDoubleDerivativeP = HPDoubleDerivative;
    }

    return HPDoubleDerivative;

}

double ManyElectrons::computeHermitePolynomialAlphaDerivative(int nValue, double position) {
    // Computes Hermite polynomials differentiated w.r.t. alpha.
    double alpha = m_parameters[0];
    double alphaSqrt = sqrt(alpha);
    double omegaSqrt = sqrt(m_omega);

    double HPDerivativePP = 0;              // d/dα H_{n-2}
    double HPDerivativeP = 0;               // d/dα H_{n-1}
    double HPDerivative = HPDerivativeP;    // d/dα H_n

    for (int n=1; n <= nValue; n++) {
        HPDerivative = omegaSqrt/alphaSqrt*position*computeHermitePolynomial(n-1, position)
                      +2*alphaSqrt*omegaSqrt*position*HPDerivativeP
                      -2*(n-1)*HPDerivativePP;
        HPDerivativePP = HPDerivativeP;
        HPDerivativeP = HPDerivative;
    }

    return HPDerivative;
}

void ManyElectrons::setUpSlaterDet() {
    // Function for setting up the Slater determinant at the begining of the simulation.
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

    m_a = zeros<mat>(m_numberOfParticles, m_numberOfParticles);
    int half = m_halfNumberOfParticles;
    for (int i=0; i < m_numberOfParticles; i++) {
        for (int j=0; j < m_numberOfParticles; j++) {
            if ( ((i < half) && (j < half)) || ((i >= half) && (j >= half)) ) { m_a(i,j) = 1./3; }
            else { m_a(i,j) = 1.; }
        }
    }

    m_spinUpSlater = zeros<mat>(m_halfNumberOfParticles, m_halfNumberOfParticles);
    m_spinDownSlater = zeros<mat>(m_halfNumberOfParticles, m_halfNumberOfParticles);

    for (int i=0; i < m_halfNumberOfParticles; i++) {
        for (int j=0; j < m_halfNumberOfParticles; j++) {
            nx = m_quantumNumbers(j, 0);
            ny = m_quantumNumbers(j, 1);
            double xSpinUp = m_system->getInitialState()->getParticles()[i]->getPosition()[0];
            double ySpinUp = m_system->getInitialState()->getParticles()[i]->getPosition()[1];
            double xSpinDown = m_system->getInitialState()->getParticles()[i+m_halfNumberOfParticles]->getPosition()[0];
            double ySpinDown = m_system->getInitialState()->getParticles()[i+m_halfNumberOfParticles]->getPosition()[1];
            m_spinUpSlater(i,j) = evaluateSingleParticleWF(nx, ny, xSpinUp, ySpinUp);
            m_spinDownSlater(i,j) = evaluateSingleParticleWF(nx, ny, xSpinDown, ySpinDown);
        }
    }

    m_spinUpSlaterInverse = m_spinUpSlater.i();
    m_spinDownSlaterInverse = m_spinDownSlater.i();
}

void ManyElectrons::updateSlaterDet(int randomParticle) {
    // Function for updating the Slater determinant after every accepted metropolis step.
    int i = randomParticle;
    std::vector<double> r_i = m_system->getParticles()[i]->getPosition();

    if (i < m_halfNumberOfParticles) {
        mat spinUpSlaterInverseOld = m_spinUpSlaterInverse;
        for (int j=0; j < i; j++) {
            double sum = 0;
            for (int l=0; l <m_halfNumberOfParticles; l++) {
                int nx = m_quantumNumbers(l, 0);
                int ny = m_quantumNumbers(l, 1);
                sum += evaluateSingleParticleWF(nx, ny, r_i[0], r_i[1])
                      *spinUpSlaterInverseOld(l,j);
            }
            for (int k=0; k < m_halfNumberOfParticles; k++) {
                m_spinUpSlaterInverse(k,j) = spinUpSlaterInverseOld(k,j)
                                            -(sum/m_ratioSlaterDet)*spinUpSlaterInverseOld(k,i);
            }
        }
        for (int j=i+1; j < m_halfNumberOfParticles; j++) {
            double sum = 0;
            for (int l=0; l <m_halfNumberOfParticles; l++) {
                int nx = m_quantumNumbers(l, 0);
                int ny = m_quantumNumbers(l, 1);
                sum += evaluateSingleParticleWF(nx, ny, r_i[0], r_i[1])
                      *spinUpSlaterInverseOld(l,j);
            }
            for (int k=0; k < m_halfNumberOfParticles; k++) {
                m_spinUpSlaterInverse(k,j) = spinUpSlaterInverseOld(k,j)
                                            -(sum/m_ratioSlaterDet)*spinUpSlaterInverseOld(k,i);
            }
        }
        for (int k=0; k < m_halfNumberOfParticles; k++) {
            m_spinUpSlaterInverse(k,i) = spinUpSlaterInverseOld(k,i)/m_ratioSlaterDet;
        }
    }
    else {
        double iHalf = i-m_halfNumberOfParticles;
        mat spinDownSlaterInverseOld = m_spinDownSlaterInverse;
        for (int j=0; j < iHalf; j++) {
            double sum = 0;
            for (int l=0; l < m_halfNumberOfParticles; l++) {
                int nx = m_quantumNumbers(l, 0);
                int ny = m_quantumNumbers(l, 1);
                sum += evaluateSingleParticleWF(nx, ny, r_i[0], r_i[1])
                      *spinDownSlaterInverseOld(l,j);
            }
            for (int k=0; k < m_halfNumberOfParticles; k++) {
                m_spinDownSlaterInverse(k,j) = spinDownSlaterInverseOld(k,j)
                                              -(sum/m_ratioSlaterDet)
                                               *spinDownSlaterInverseOld(k, iHalf);
            }
        }
        for (int j=iHalf+1; j < m_halfNumberOfParticles; j++) {
            double sum = 0;
            for (int l=0; l < m_halfNumberOfParticles; l++) {
                int nx = m_quantumNumbers(l, 0);
                int ny = m_quantumNumbers(l, 1);
                sum += evaluateSingleParticleWF(nx, ny, r_i[0], r_i[1])
                      *spinDownSlaterInverseOld(l,j);
            }
            for (int k=0; k < m_halfNumberOfParticles; k++) {
                m_spinDownSlaterInverse(k,j) = spinDownSlaterInverseOld(k,j)
                                              -(sum/m_ratioSlaterDet)
                                               *spinDownSlaterInverseOld(k, iHalf);
            }
        }
        for (int k=0; k < m_halfNumberOfParticles; k++) {
            m_spinDownSlaterInverse(k, iHalf) = spinDownSlaterInverseOld(k, iHalf)/m_ratioSlaterDet;
        }
    }
}
