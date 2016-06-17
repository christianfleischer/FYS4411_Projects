#ifndef PROJECT2_HAMILTONIAN_H
#define PROJECT2_HAMILTONIAN_H
#include <vector>

class Hamiltonian {
public:
    Hamiltonian(class System* system, bool analyticalKinetic);
    double computeKineticEnergy(std::vector<class Particle*> particles);
    virtual std::vector<double> computeLocalEnergy(std::vector<class Particle*> particles) = 0;
    double getAnalytic(){ return m_analyticalKinetic; }

protected:
    class System* m_system = nullptr;
    bool m_analyticalKinetic = false;
};

#endif // PROJECT2_HAMILTONIAN_H
