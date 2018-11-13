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
    std::string fileName = "16O_8_8_1.30.txt";
    std::string filePotential = "KS-POTENTIAL.dat";
    std::string fileDensity = "208Pb(p)_SKP_siCoul_density-2Rn.txt";
    Parameters::ElementConstants::initialiseElementConstants(fileName);
    const double H = 0.1;
    const double xin = 1e-6;
    const double xfin = 3*Parameters::ElementConstants::Rn();
    int potType = 1;	// set: 0: Total potential | 1: Woods-Saxon | 2: Harmonic Oscillator | 3: Test potential | 4: KS/Other potential
    int empyDensType = 0;	// set: 0: Protons SoG | 1: Neutrons SoG | 2: MC/Other density
    Parameters::IntegrationParameters::initialiseIntegrationParameters(xin, xfin);
    Parameters::NucleonType::initialiseNucleonType(true); 	// false: run simulation for protons | true: run simulation for neutrons
    Parameters::SpinOrbit::initialiseSpinOrbit(false); 		// false: Spin-Orbit OFF | true: Spin-Orbit ON
    
    mkdir("Outputs", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH); // Create outputs folder
    mkdir("Outputs/PotStep", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH); // Create potential step folder
    std::string inputPath = "Inputs/";
    std::string outputPath = "Outputs/";
    std::string stepPath = "Outputs/PotStep/";
    clock_t start = clock(); // Start time

    std::string fileQuantNumb;
    int nMatRows;
    if (!Parameters::SpinOrbit::SpinOrbitOn())
    {
    	fileQuantNumb = "orbitals.txt";
    	nMatRows = 36;
    }
    else
    {
    	fileQuantNumb = "orbitals-so.txt";
    	nMatRows = 28;
    }

    //Load the matrix of quantum numbers for each state from "orbitals.txt" and SoG parameters from file
    std::vector <std::vector <double>> orbitals (nMatRows);
    std::ifstream in (inputPath + fileQuantNumb);
    readMatrix(orbitals, in, false);
    in.close();

    in.open(inputPath + fileName);
    std::vector<std::vector<double>> qrParam(12);
    readMatrix(qrParam, in, false);
    in.close();

    Element nuclei(orbitals);

	double nclMass = Parameters::mp;
	if(Parameters::NucleonType::isNeutron())
		nclMass = Parameters::mn;

	//Test theoretical and sog density for initial potential
	std::unique_ptr<InitialPot> pot_ptr;
	switch(potType)
	{
		case 0:
			//pot_ptr = std::make_unique<TotPot> (nclMass);
		break;
		case 1:
			pot_ptr = std::make_unique<WSaxPot> (Parameters::ElementConstants::Rn(), Parameters::a0, nclMass);
		break;
		case 2:
			pot_ptr = std::make_unique<HOPot> (nclMass);
		break;
		case 3:
			pot_ptr = std::make_unique<TestPot> (nclMass);
		break;
		case 4:

			in.open(inputPath + filePotential);
		    std::map<double, double> initialPotential;
		    readMap(initialPotential, in, false);

			pot_ptr = std::make_unique<PotOut> (initialPotential, nclMass);
		break;
	}

	//std::string potTot = "potTot" + potType;
    std::map<double, double> inPotential;
    const unsigned long NSteps = std::abs(Parameters::IntegrationParameters::x1() - Parameters::IntegrationParameters::x0())/H;
    double rad = Parameters::IntegrationParameters::x0();
    for (int i = 0; i < NSteps +1; ++i)
    {
 		inPotential[rad] = pot_ptr -> potential(rad);
 		rad += H;
    }
    std::ofstream fOut(outputPath + "refInitialPotential.txt");
    writeMap(inPotential, fOut, false);
    fOut.close();


    const Schroddy sh(*pot_ptr, H);
    ElementEigenfunctions elEigf = nuclei.orbitalEigenfunction(sh, orbitals);

    std::unique_ptr<NuclearDensity> NDens_ptr;
    switch(empyDensType)
    {
    	case 0:
    		NDens_ptr = std::make_unique<NuclearDensityWithSOG>();
    		NDens_ptr -> theoreticalDensity(elEigf, nuclei.getLevelDegeneration());
    		NDens_ptr -> benchmarkDensity(qrParam, H);

    	break;
    	case 1:
    		NDens_ptr = std::make_unique<NuclearDensityWithNeutronsSOG>();
    		NDens_ptr -> theoreticalDensity(elEigf, nuclei.getLevelDegeneration());
    		NDens_ptr -> benchmarkDensity(qrParam, H);

    	break;
    	case 2:
    		NDens_ptr = std::make_unique<NuclearDensityWithMC>();
    		NDens_ptr -> theoreticalDensity(elEigf, nuclei.getLevelDegeneration());
    		in.open(inputPath + fileDensity);
    		std::vector<std::vector<double>> mcDensity(NDens_ptr -> getTheoreticalDensity().size());
    		readMatrix(mcDensity, in, false);
    		in.close();
    		NDens_ptr -> benchmarkDensity(mcDensity);

    	break;
    }

    fOut.open(outputPath + "refInitialDensity.txt");
    fOut << NDens_ptr -> getTheoreticalDensity();
    fOut.close();
    fOut.open(outputPath + "refSogDensity.txt");
    fOut << NDens_ptr -> getBenchmarkDensity();
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
    KohnShamInverseWithLB ksi(*pot_ptr, H);
    KohnShamInverseWithLB tempKsi = ksi;
    /*KohnShamInverseWithJW ksi(potTot, H);
    KohnShamInverseWithJW tempKsi = ksi;*/
    /*KohnShamInverseWithWP ksi(potTot, H, nuclei.getLevelEigenvalue(), elEigf);
    KohnShamInverseWithWP tempKsi = ksi;*/
    ksi.KSinverse(*NDens_ptr, tempKsi);
    fOut << ksi.getKSPot();
    fOut.close();

    unsigned long loops = 0;
    while (!NDens_ptr -> hasConverged()) //simulation loop
    {
        const PotOut po(ksi, Parameters::mn, 0);
        const Schroddy sh_(po, H);
        elEigf = nuclei.orbitalEigenfunction(sh_, orbitals);
        NDens_ptr -> theoreticalDensity(elEigf, nuclei.getLevelDegeneration());
        KohnShamInverseWithLB tempKsi_ = ksi;
        //KohnShamInverseWithJW tempKsi_ = ksi;
        /*ksi.setElementEigenvalues(nuclei.getLevelEigenvalue());
        ksi.setElementEigenfunctions(elEigf);
        KohnShamInverseWithWP tempKsi_ = ksi;*/
        ksi.KSinverse(*NDens_ptr, tempKsi_);

        std::string nfile = std::to_string(loops);
        KSPotential finalPotential = ksi.getKSPot();
        fOut.open(stepPath + nfile + "StepPotential.txt");
        fOut << finalPotential;
        fOut.close();

        ++loops;
        //if (loops%10 == 0)
            std::cerr << "Convergence distance: " << NDens_ptr -> distanceToConvergence() << " with epsilon "
            << NDens_ptr -> epsilon() << " after " << loops << " iterations. " << std::endl;
    }

    EnergyTOT totalenergy;
    fOut.open(outputPath + "refTotalEnergy.txt");
    fOut << totalenergy.energyTot(elEigf, nuclei.getLevelEigenvalue());
    fOut.close();
    fOut.open(outputPath + "refFinalDensity.txt");
    fOut << NDens_ptr -> getTheoreticalDensity();
    fOut.close();
    //NDens.densityError();
    fOut.open(outputPath + "refDensityError.txt");
    fOut << NDens_ptr -> densityError();
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

    RecapFile rcfile(H, loops, seconds, fileName, filePotential, fileDensity, fileQuantNumb, potType, empyDensType);
    rcfile.paramRecap();

    std::cout << "CONVERGENCE IS ACHIEVED in: " << seconds << " seconds!" << " GREAT JOB!" << std::endl;
    std::cout << "Program executed successfully." << std::endl;
    return 0;
}
