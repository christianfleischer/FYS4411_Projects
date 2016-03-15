#pragma once

class SteepestDescent {
public:
    SteepestDescent(class System* system, double stepLengthSD);
    void obtainOptimalAlpha(double alpha, double tol, int maxIterations,
                            int numberOfMetropolisSteps, bool importanceSampling);

private:
    double m_stepLengthSD = 0;
    class System* m_system = nullptr;
};
