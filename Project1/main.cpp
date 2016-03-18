#include <iostream>
#include "system.h"
#include "sampler.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "WaveFunctions/simplegaussian.h"
#include "WaveFunctions/repulsivegaussian.h"
#include "Hamiltonians/hamiltonian.h"
#include "Hamiltonians/harmonicoscillator.h"
#include "Hamiltonians/harmonicoscillatorrepulsive.h"
#include "InitialStates/initialstate.h"
#include "InitialStates/randomuniform.h"
#include "VariationMethods/steepestdescent.h"
#include "Math/random.h"
#include <mpi.h>

using namespace std;


int main(/*int nargs, char* args[]*/) {

//    int numprocs, my_rank;
//    MPI_Init(&nargs, &args);
//    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
//    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    int numberOfDimensions  = 3;
    int numberOfParticles   = 10;
    int numberOfSteps       = (int) 1e6;    // Monte Carlo cycles
    double omega            = 1.0;          // Oscillator frequency.
    double alpha            = 0.5;          // Variational parameter.
    double beta             = 2.82843;      // Variational parameter.
    double gamma            = 2.82843;
    double a                = 0.0043;       // Hard core boson diameter.
    double stepLength       = 0.1;          // Metropolis step length.
    double equilibration    = 0.1;          // Amount of the total steps used for equilibration.
    double dt               = 0.01;         // Time step for importance sampling.
    bool analyticalKinetic  = true;
    bool importanceSampling = false;
    bool repulsion          = true;         // Switch for interacting system or not.
    bool saveEnergies       = false;
    bool savePositions      = false;
    bool showProgress       = true;
    bool printToTerminal    = true;

    int num_my_cycles = numberOfSteps;// /numprocs;

//    cout << "  -- Settings -- " << boolalpha << endl;
//    cout << " Analytical Kinetic : " << analyticalKinetic << endl;
//    cout << " Importance Sampling : " << importanceSampling << endl;
//    cout << " Repulsion : " << repulsion << endl;

    // Initiate System
    System* system = new System();
    // Select which Hamiltonian and trial wave function to use (interacting or non-interacting)
    if (repulsion) {
    system->setHamiltonian              (new HarmonicOscillatorRepulsive(system, omega, a, gamma, analyticalKinetic));
    system->setWaveFunction             (new RepulsiveGaussian(system, alpha, beta, a));
    }
    else{
    system->setHamiltonian              (new HarmonicOscillator(system, omega, analyticalKinetic));
    system->setWaveFunction             (new SimpleGaussian(system, alpha));
    }
    // RandomUniform creates a random initial state
    system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles/*, my_rank*/));
    system->setEquilibrationFraction    (equilibration);
    system->setStepLength               (stepLength);
    system->setTimeStep                 (dt);
    // Start Monte Carlo simulation
    system->runMetropolisSteps          (num_my_cycles, importanceSampling, saveEnergies,
                                         savePositions, showProgress, printToTerminal);

//    double e = system->getSampler()->getEnergy();
//    double total_e = 0;
//    MPI_Reduce(&e, &total_e, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
//    if (my_rank == 0){
//        total_e /= numprocs;
//        cout << total_e << endl;
//    }
//    MPI_Finalize();

    // Steepest descent:
    cout << "Optimizing alpha using steepest descent:" << endl;
    int maxIterations           = 1000;
    int numberOfStepsSD         = (int) 1e5;
    double stepLengthSD         = 0.01;
    double initialAlpha         = 0.7;
    double tol                  = 1e-6;//0.001;
    bool importanceSamplingSD   = false;
    SteepestDescent* steepestDescent = new SteepestDescent(system, stepLengthSD);
    steepestDescent->obtainOptimalAlpha(initialAlpha, tol, maxIterations, numberOfStepsSD, importanceSamplingSD);


    return 0;
}

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
