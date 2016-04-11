#pragma once
#include <vector>

class Hamiltonian {
public:
    Hamiltonian(class System* system, bool analyticalKinetic);
    double computeKineticEnergy(std::vector<class Particle*> particles);
    virtual double computeLocalEnergy(std::vector<class Particle*> particles) = 0;
    double getAnalytic(){ return m_analyticalKinetic; }

protected:
    class System* m_system = nullptr;
    bool m_analyticalKinetic = false;
};

