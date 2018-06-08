#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "includes.h"


int main(int argc, const char * argv[])
{
    std::string fileName = "40Ca_20_20_1.45.txt";
    std::string filePotential = "PROTONI_208pb_KS-POTENTIAL.dat";
    Parameters::ElementConstants::initialiseElementConstants(fileName);
    const double H = 0.1;
    const double xin = 1e-6;
    const double xfin = 3*Parameters::ElementConstants::Rn();
    Parameters::IntegrationParameters::initialiseIntegrationParameters(xin, xfin);
    Parameters::NucleonType::initialiseNucleonType(false); 	//false: run simulation for protons
                                                            //true: run simulation for neutrons
    
    mkdir("Outputs", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH); // Create outputs folder
    std::string inputPath = "Inputs/";
    std::string outputPath = "Outputs/";
    clock_t start = clock(); // Start time

    //Load the matrix of quantum numbers for each state from "orbitals.txt"
    std::vector <std::vector <double>> orbitals (36);
    std::ifstream in (inputPath + "orbitals.txt");
    readMatrix(orbitals, in, false);
    in.close();

    //in.open(inputPath + "16O.txt");
    in.open(inputPath + fileName);
    std::vector<std::vector<double>> qrParam(12);
    readMatrix(qrParam, in, false);
    in.close();

    //load initial potential from file
    /*in.open(inputPath + filePotential);
    std::map<double, double> initialPotential;
    readMap(initialPotential, in, false);*/

    Element nuclei(orbitals);

	double nclMass = Parameters::mp;
	if(Parameters::NucleonType::isNeutron())
		nclMass = Parameters::mn;

    //Test theoretical and sog density for initial harmonic potential
    //const TotPot potTot (nclMass);
    const WSaxPot potTot (Parameters::ElementConstants::Rn(), Parameters::a0, nclMass);
    //const HOPot potTot(nclMass); //default is ground state
    //const TestPot potTot(nclMass);
    //const PotOut potTot (initialPotential, nclMass);

    std::map<double, double> inPotential;
    const unsigned long NSteps = std::abs(Parameters::IntegrationParameters::x1() - Parameters::IntegrationParameters::x0())/H;
    double rad = Parameters::IntegrationParameters::x0();
    for (int i = 0; i < NSteps +1; ++i)
    {
 		inPotential[rad] = potTot.potential(rad);
 		rad += H;
    }
    std::ofstream fOut(outputPath + "refInitialPotential.txt");
    writeMap(inPotential, fOut, false);
    fOut.close();


    const Schroddy sh(potTot, H);
    ElementEigenfunctions elEigf = nuclei.orbitalEigenfunction(sh, orbitals);
    NuclearDensityWithSOG NDens;
    //NuclearDensityWithNeutronsSOG NDens;
    //NuclearDensityWithMC NDens;
    NDens.theoreticalDensity(elEigf, nuclei.getLevelDegeneration());
    NDens.benchmarkDensity(qrParam, H);

    /*in.open(inputPath + "HODensity-20NN-3Rn.txt");
    std::vector<std::vector<double>> mcDensity(NDens.getTheoreticalDensity().size());
    readMatrix(mcDensity, in, false);
    in.close();
    NDens.benchmarkDensity(mcDensity);*/

    fOut.open(outputPath + "refInitialDensity.txt");
    fOut << NDens.getTheoreticalDensity();
    fOut.close();
    fOut.open(outputPath + "refSogDensity.txt");
    fOut << NDens.getBenchmarkDensity();
    fOut.close();
    ElementEigenvalues initialEigenvalues = nuclei.getLevelEigenvalue();
    fOut.open(outputPath + "refInitialEigenvalues.txt");
    writeMatrix(initialEigenvalues, fOut, false);
    fOut.close();
    fOut.open(outputPath + "refInitialEigenfunctions.txt");
    writeElementEigenfunctions(elEigf, fOut);
    fOut.close();
    /*EnergyTOT totalenergy;
    fOut.open(outputPath + "refTotalEnergy.txt");
    fOut << totalenergy.energyTot(elEigf, nuclei.getLevelEigenvalue(), H);
    fOut.close();*/

    //Test Kohn-Sham inversion for initial harmonic potential
    fOut.open(outputPath + "refFirstKSPotential.txt");
    KohnShamInverseWithLB ksi(potTot, H);
    KohnShamInverseWithLB tempKsi = ksi;
    /*KohnShamInverseWithJW ksi(potTot, H);
    KohnShamInverseWithJW tempKsi = ksi;*/
    /*KohnShamInverseWithWP ksi(potTot, H, nuclei.getLevelEigenvalue(), elEigf);
    KohnShamInverseWithWP tempKsi = ksi;*/
    ksi.KSinverse(NDens, tempKsi);
    fOut << ksi.getKSPot();
    fOut.close();

    unsigned long loops = 0;
    while (!NDens.hasConverged()) //simulation loop
    {
        const PotOut po(ksi, Parameters::mn, 0);
        const Schroddy sh_(po, H);
        elEigf = nuclei.orbitalEigenfunction(sh_, orbitals);
        NDens.theoreticalDensity(elEigf, nuclei.getLevelDegeneration());
        KohnShamInverseWithLB tempKsi_ = ksi;
        //KohnShamInverseWithJW tempKsi_ = ksi;
        /*ksi.setElementEigenvalues(nuclei.getLevelEigenvalue());
        ksi.setElementEigenfunctions(elEigf);
        KohnShamInverseWithWP tempKsi_ = ksi;*/
        ksi.KSinverse(NDens, tempKsi_);

        ++loops;
        //if (loops%10 == 0)
            std::cerr << "Convergence distance: " << NDens.distanceToConvergence() << " with epsilon "
            << NDens.epsilon() << " after " << loops << " iterations. " << std::endl;
    }

    EnergyTOT totalenergy;
    fOut.open(outputPath + "refTotalEnergy.txt");
    fOut << totalenergy.energyTot(elEigf, nuclei.getLevelEigenvalue());
    fOut.close();
    fOut.open(outputPath + "refFinalDensity.txt");
    fOut << NDens.getTheoreticalDensity();
    fOut.close();
    //NDens.densityError();
    fOut.open(outputPath + "refDensityError.txt");
    fOut << NDens.densityError();
    fOut.close();
    ElementEigenvalues finalEigenvalues = nuclei.getLevelEigenvalue();
    fOut.open(outputPath + "refFinalEigenvalues.txt");
    writeMatrix(finalEigenvalues, fOut, false);
    fOut.close();
    fOut.open(outputPath + "refFinalEigenfunctions.txt");
    writeElementEigenfunctions(elEigf, fOut);
    fOut.close();
    KSPotential finalPotential = ksi.getKSPot();
    fOut.open(outputPath + "refFinalPotential.txt");
    fOut << finalPotential;
    fOut.close();

    clock_t end = clock(); // Finish time
    double seconds = (((double)(end - start))/CLOCKS_PER_SEC);

    std::cout << "CONVERGENCE IS DONE in: " << seconds << " seconds!" << " GREAT JOB!" << std::endl;
    std::cout << "Program executed successfully." << std::endl;
    return 0;
}
