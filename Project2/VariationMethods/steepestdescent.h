#pragma once
#include "iostream"

class SteepestDescent {
public:
    SteepestDescent(class System* system, double stepLengthSD);
    void obtainOptimalParameter(double parameter, std::string parameterName, double tol,
                                int maxIterations, int numberOfMetropolisSteps, bool importanceSampling);

private:
    double m_stepLengthSD = 0;
    class System* m_system = nullptr;
};
