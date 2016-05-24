#ifndef PROJECT2_WAVEFUNCTION_H
#define PROJECT2_WAVEFUNCTION_H
#include <vector>


class WaveFunction {
public:
    WaveFunction(class System* system);
    int     getNumberOfParameters() { return m_numberOfParameters; }
    std::vector<double> getParameters() { return m_parameters; }
    void adjustParameter(double parameter, int parameterNumber);
    virtual double evaluate(std::vector<class Particle*> particles) = 0;
    virtual double computeDoubleDerivative(std::vector<class Particle*> particles) = 0;
    virtual double computeMetropolisRatio(std::vector<class Particle*> particles, int randomParticle,
                                          std::vector<double> positionChange) = 0;
    virtual std::vector<double> computeDerivative(std::vector<class Particle*> particles,
                                                  int randomParticle) = 0;
    virtual std::vector<double> computeDerivativeWrtParameters(std::vector<Particle *> particles) = 0;
    virtual void updateSlaterDet(int randomParticle) = 0;

protected:
    int     m_numberOfParameters = 0;
    std::vector<double> m_parameters = std::vector<double>();
    class System* m_system = nullptr;
};

#endif // PROJECT2_WAVEFUNCTION_H
