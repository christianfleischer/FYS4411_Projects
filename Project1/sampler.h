#pragma once

class Sampler {
public:
    Sampler(class System* system);
    void setNumberOfMetropolisSteps(int steps);
    void sample(bool acceptedStep, bool saveEnergies, bool savePositions);
    void printOutputToTerminal();
    void computeAverages();
    void saveToFile(double localEnergy, bool saveEnergies, bool savePositions);
    double getEnergy()              { return m_energy; }
    double getWaveFuncDerivative()  { return m_waveFuncDerivative; }
    double getWaveFuncEnergy()      { return m_waveFuncEnergy; }

private:
    int     m_numberOfMetropolisSteps = 0;
    int     m_stepNumber = 0;
    int     m_cumulativeAcceptedSteps = 0;
    double  m_energy = 0;
    double  m_squaredEnergy = 0;
    double  m_cumulativeEnergy = 0;
    double  m_cumulativeSquaredEnergy = 0;
    double  m_variance = 0;
    double  m_acceptanceRate = 0;
    double  m_waveFuncDerivative = 0;
    double  m_waveFuncEnergy = 0;
    double  m_cumulativeWFuncDerivative = 0;
    double  m_cumulativeWFuncEnergy = 0;
    class System* m_system = nullptr;
};
