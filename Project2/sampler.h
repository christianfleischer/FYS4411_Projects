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
    double getWaveFuncDerivativeAlpha()  { return m_waveFuncDerivativeAlpha; }
    double getWaveFuncEnergyAlpha()      { return m_waveFuncEnergyAlpha; }
    double getWaveFuncDerivativeBeta()  { return m_waveFuncDerivativeBeta; }
    double getWaveFuncEnergyBeta()      { return m_waveFuncEnergyBeta; }

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
    double  m_waveFuncDerivativeAlpha = 0;
    double  m_waveFuncEnergyAlpha = 0;
    double  m_cumulativeWFuncDerivativeAlpha = 0;
    double  m_cumulativeWFuncEnergyAlpha = 0;
    double  m_waveFuncDerivativeBeta = 0;
    double  m_waveFuncEnergyBeta = 0;
    double  m_cumulativeWFuncDerivativeBeta = 0;
    double  m_cumulativeWFuncEnergyBeta = 0;
    class System* m_system = nullptr;
};
