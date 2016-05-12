#include <iostream>
#include <fstream>
#include "system.h"
#include "sampler.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "WaveFunctions/simplegaussian.h"
#include "WaveFunctions/repulsivegaussian.h"
#include "WaveFunctions/twoelectrons.h"
#include "WaveFunctions/manyelectrons.h"
#include "Hamiltonians/hamiltonian.h"
#include "Hamiltonians/harmonicoscillator.h"
#include "Hamiltonians/harmonicoscillatorrepulsive.h"
#include "Hamiltonians/harmonicoscillatorelectrons.h"
#include "InitialStates/initialstate.h"
#include "InitialStates/randomuniform.h"
#include "VariationMethods/steepestdescent.h"
#include "Math/random.h"
#include <mpi.h>

using namespace std;


int main(int nargs, char* args[]) {

    int numprocs, my_rank;
    double totalE, totalVariance, totalAcceptanceRate, finalMeanDistance;
    double timeStart, timeEnd, totalTime;

    MPI_Init(&nargs, &args);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    timeStart = MPI_Wtime();

    int numberOfDimensions  = 2;
    int numberOfParticles   = 2;
    int numberOfSteps       = (int) 1e5;    // Monte Carlo cycles
    double omega            = 1.0;          // Oscillator frequency.
    double alpha            = 0.737505;//0.7;          // Variational parameter.
    double beta             = 0.364;//0.506577;//2.82843;      // Variational parameter.
    double gamma            = 2.82843;
    double a                = 0.0043;       // Hard core boson diameter.
    double stepLength       = 0.1;          // Metropolis step length.
    double equilibration    = 0.1;          // Amount of the total steps used for equilibration.
    double dt               = 0.01;         // Time step for importance sampling.
    double aElectrons       = 1; //1./3
    double C                = 1;            // Norm constant.
    bool analyticalKinetic  = true;
    bool importanceSampling = true;
    bool repulsion          = false;         // Switch for interacting system or not.
    bool saveEnergies       = false;
    bool savePositions      = false;
    bool showProgress       = true;
    bool printToTerminal    = true;
    bool quantumDots        = true;

    int numMyCycles = numberOfSteps/numprocs;

//    cout << "  -- Settings -- " << boolalpha << endl;
//    cout << " Analytical Kinetic : " << analyticalKinetic << endl;
//    cout << " Importance Sampling : " << importanceSampling << endl;
//    cout << " Repulsion : " << repulsion << endl;

    // Initiate System
    System* system = new System();
    // RandomUniform creates a random initial state
    system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles, my_rank));
    // Select which Hamiltonian and trial wave function to use (interacting or non-interacting)
    if (repulsion && !quantumDots) {
    system->setHamiltonian              (new HarmonicOscillatorRepulsive(system, omega, a, gamma, analyticalKinetic));
    system->setWaveFunction             (new RepulsiveGaussian(system, alpha, beta, a));
    }
    if (!repulsion && !quantumDots) {
    system->setHamiltonian              (new HarmonicOscillator(system, omega, analyticalKinetic));
    system->setWaveFunction             (new SimpleGaussian(system, alpha));
    }
    if (quantumDots) {
        if (numberOfParticles == 1) {
            system->setHamiltonian      (new HarmonicOscillatorElectrons(system, omega, analyticalKinetic, repulsion));
            system->setWaveFunction     (new TwoElectrons(system, alpha, beta, omega, aElectrons, C));
        }
        else {
            system->setHamiltonian      (new HarmonicOscillatorElectrons(system, omega, analyticalKinetic, repulsion));
            system->setWaveFunction     (new ManyElectrons(system, alpha, beta, omega, C));
        }
    }
    system->setEquilibrationFraction    (equilibration);
    system->setStepLength               (stepLength);
    system->setTimeStep                 (dt);
    system->setMyRank                   (my_rank);
    // Optimize parameters
    //system->optimizeParameters          (system, alpha, beta);
    system->setSaveEnergies             (saveEnergies);
    system->setSavePositions            (savePositions);
    // Start Monte Carlo simulation
    system->runMetropolisSteps          (numMyCycles, importanceSampling, showProgress, printToTerminal);
    // Compute MPI averages etc.
    system->MPI_CleanUp                 (totalE, totalVariance, totalAcceptanceRate, finalMeanDistance,
                                         timeStart, timeEnd, totalTime, numprocs, numberOfSteps);
    // Merge the files from the nodes into one data file
    system->mergeOutputFiles            (numprocs);

    //ManyElectrons* manyelectrons = new ManyElectrons(system, 1, 1, 1, 1, 1);
    //cout << manyelectrons->computeHermitePolynomialDoubleDerivative(4, 2) << endl;

    return 0;
}

//12: 65.7
//20: 155.868

/*
 -- System info --
 Number of particles  : 50
 Number of dimensions : 3
 Number of Metropolis steps run : 10^6
 Number of equilibration steps  : 10^5

  -- Wave function parameters --
 Number of parameters : 2
 Parameter 1 : 0.5
 Parameter 2 : 2.82843

  -- Results --
 Energy : 127.306
 Variance : 0.191128
 Acceptance Rate : 0.891037

Computation Time : 13577

Optimizing alpha using steepest descent:

 -- System info --
 Number of particles  : 500
 Number of dimensions : 3
 Number of Metropolis steps run : 10^5
 Number of equilibration steps  : 10^4

  -- Wave function parameters --
 Number of parameters : 1
 Parameter 1 : 0.5

  -- Results --
 Energy : 750
 Variance : 7.64674e-06
 Acceptance Rate : 0.911399

Computation Time : 4580.39


*/
