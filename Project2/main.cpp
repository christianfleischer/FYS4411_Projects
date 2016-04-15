#include <iostream>
#include <fstream>
#include "system.h"
#include "sampler.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "WaveFunctions/simplegaussian.h"
#include "WaveFunctions/repulsivegaussian.h"
#include "WaveFunctions/twoelectrons.h"
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
    double total_e, total_variance, total_acceptanceRate, final_mean_distance;
    double time_start, time_end, total_time;

    MPI_Init(&nargs, &args);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    time_start = MPI_Wtime();

    int numberOfDimensions  = 2;
    int numberOfParticles   = 2;
    int numberOfSteps       = (int) 1e6;    // Monte Carlo cycles
    double omega            = 1.0;          // Oscillator frequency.
    double alpha            = 0.737503;          // Variational parameter.
    double beta             = 0.506579;      // Variational parameter.
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
    bool saveEnergies       = true;
    bool savePositions      = false;
    bool showProgress       = true;
    bool printToTerminal    = true;
    bool quantumDots        = true;

    int num_my_cycles = numberOfSteps/numprocs;

//    cout << "  -- Settings -- " << boolalpha << endl;
//    cout << " Analytical Kinetic : " << analyticalKinetic << endl;
//    cout << " Importance Sampling : " << importanceSampling << endl;
//    cout << " Repulsion : " << repulsion << endl;

    // Initiate System
    System* system = new System();
    // Select which Hamiltonian and trial wave function to use (interacting or non-interacting)
    if (repulsion && !quantumDots) {
    system->setHamiltonian              (new HarmonicOscillatorRepulsive(system, omega, a, gamma, analyticalKinetic));
    system->setWaveFunction             (new RepulsiveGaussian(system, alpha, beta, a));
    }
    if (!repulsion && !quantumDots){
    system->setHamiltonian              (new HarmonicOscillator(system, omega, analyticalKinetic));
    system->setWaveFunction             (new SimpleGaussian(system, alpha));
    }
    if (quantumDots){
        if (numberOfParticles == 2){
            system->setHamiltonian              (new HarmonicOscillatorElectrons(system, omega, analyticalKinetic));
            system->setWaveFunction             (new TwoElectrons(system, alpha, beta, omega, aElectrons, C));
        }
    }
    // RandomUniform creates a random initial state
    system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles, my_rank));
    system->setEquilibrationFraction    (equilibration);
    system->setStepLength               (stepLength);
    system->setTimeStep                 (dt);
    system->setMyRank                   (my_rank);
    system->setSaveEnergies             (saveEnergies);
    // Start Monte Carlo simulation
    system->runMetropolisSteps          (num_my_cycles, importanceSampling, saveEnergies,
                                         savePositions, showProgress, printToTerminal);

    double e = system->getSampler()->getEnergy();
    double variance = system->getSampler()->getVariance();
    double acceptanceRate = system->getSampler()->getAcceptanceRate();
    double mean_distance = system->getSampler()->getMeanDistance();

    MPI_Reduce(&e, &total_e, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&variance, &total_variance, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&acceptanceRate, &total_acceptanceRate, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&mean_distance, &final_mean_distance, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    time_end = MPI_Wtime();
    total_time = time_end-time_start;

    if (my_rank == 0){
        total_e /= numprocs;
        total_variance /= numprocs;
        total_acceptanceRate /= numprocs;
        final_mean_distance /= numprocs;
        system->setNumberOfMetropolisSteps(numberOfSteps);
        system->setComputationTime(total_time);
        system->getSampler()->setEnergy(total_e);
        system->getSampler()->setVariance(total_variance);
        system->getSampler()->setAcceptanceRate(total_acceptanceRate);
        system->getSampler()->setMeanDistance(final_mean_distance);
        system->getSampler()->printOutputToTerminal();
    }
    if (saveEnergies) fclose(system->getEnergiesFile());
    MPI_Finalize();

    if (saveEnergies){
        std::ofstream outfile("energies.dat", std::ios_base::binary);
        for (int i=0; i < numprocs; i++){
            char nodeFileName[100];
            sprintf(nodeFileName, "energiesNode%i.dat", i);
            std::ifstream nodeFile(nodeFileName, std::ios_base::binary);

            outfile << nodeFile.rdbuf();
            nodeFile.close();
            remove(nodeFileName);
        }
        outfile.close();
    }

    /*
    // Steepest descent:
    cout << "Optimizing alpha using steepest descent:" << endl;
    int maxIterations             = 100;
    int numberOfStepsSD           = (int) 1e5;
    double stepLengthSD           = 0.01;
    double initialAlpha           = 0.7;
    double tol                    = 1e-6;//0.001;
    bool importanceSamplingSD     = false;
    std::string parameterAlpha    = "alpha";

    SteepestDescent* steepestDescent = new SteepestDescent(system, stepLengthSD);
    steepestDescent->obtainOptimalParameter(initialAlpha, parameterAlpha, tol, maxIterations,
                                            numberOfStepsSD, importanceSamplingSD);

    cout << "Optimizing beta using steepest descent:" << endl;
    maxIterations             = 100;
    numberOfStepsSD           = (int) 1e5;
    stepLengthSD              = 0.01;
    double initialBeta        = 0.505;
    tol                       = 1e-6;//0.001;
    importanceSamplingSD      = false;
    std::string parameterBeta = "beta";

    steepestDescent = new SteepestDescent(system, stepLengthSD);
    steepestDescent->obtainOptimalParameter(initialBeta, parameterBeta, tol, maxIterations,
                                            numberOfStepsSD, importanceSamplingSD);
    */

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
