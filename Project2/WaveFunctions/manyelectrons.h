#ifndef PROJECT2_MANYELECTRONS_H
#define PROJECT2_MANYELECTRONS_H
#include "wavefunction.h"
#include <armadillo>
using namespace arma;

class ManyElectrons : public WaveFunction {
public:
    ManyElectrons(class System* system, double alpha, double beta, double omega, double C);
    double evaluate(std::vector<class Particle*> particles);
    double computeDoubleDerivative(std::vector<class Particle*> particles);
    double computeDerivativeWrtAlpha(std::vector<Particle *> particles);
    double computeDerivativeWrtBeta(std::vector<Particle *> particles);
    double computeMetropolisRatio(std::vector<Particle *> particles, int randomParticle,
                                  std::vector<double> positionChange);
    double computeHermitePolynomial(int nValue, double position);
    double computeHermitePolynomialDerivative(int nValue, double position);
    double computeHermitePolynomialDoubleDerivative(int nValue, double position);
    double evaluateSingleParticleWF(int nx, int ny, double x, double y);
    double computeSPWFDoubleDerivative(int nx, int ny, double x, double y);
    std::vector<double> computeDerivative(std::vector<class Particle*> particles,
                                          int randomParticle);
    std::vector<double> computeSlaterGradient(std::vector<Particle *> particles, int i);
    std::vector<double> computeJastrowGradient(std::vector<Particle *> particles, int i);
    std::vector<double> computeSPWFDerivative(int nx, int ny, double x, double y);
    void setUpSlaterDet();
    void updateSlaterDet(int randomParticle);

private:
    int m_numberOfParticles = 0;
    int m_halfNumberOfParticles = 0;
    int m_k = 0;
    double m_omega = 0;
    double m_C = 0;
    double m_metropolisRatio = 0;
    mat m_quantumNumbers;
    mat m_spinUpSlater;
    mat m_spinDownSlater;
    mat m_spinUpSlaterInverse;
    mat m_spinDownSlaterInverse;
    mat m_a;
};

#endif // PROJECT2_MANYELECTRONS_H
