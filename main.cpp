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

    const double H = 0.1;
    mkdir("Outputs", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH); // Create outputs folder
    std::string inputPath = "Inputs/";
    std::string outputPath = "Outputs/";
    clock_t start = clock(); // Start time

    //Load the matrix of quantum numbers for each state from "orbitals.txt"
    std::vector <std::vector <unsigned int>> orbitals (36);
    std::fstream in (inputPath + "orbitals.txt", std::ios::in);
    readMatrix(orbitals, in, false);
    in.close();

    in.open(inputPath + "40Ca.txt");
    std::vector<std::vector<double>> qrParam(12);
    readMatrix(qrParam, in, false);

    Element nuclei(orbitals);

    //Test theoretical and sog density for initial harmonic potential
    const HOPot potTot(Parameters::mn); //default is ground state

    std::map<double, double> inPotential;
    const unsigned long NSteps = std::abs(Parameters::x_fin - Parameters::x_in)/H;
    double rad = Parameters::x_in;
    for (int i = 0; i < NSteps +1; ++i)
    {
 		inPotential[rad] = potTot.potential(rad);
 		rad += H;
    }
    std::ofstream fOut(outputPath + "refInitialPotential.txt");
    for (auto i: inPotential)
    	fOut << i.first << "\t" << i.second << std::endl;
    fOut.close();

    //const TotPot potTot (Parameters::mn);
    //const WSaxPot potTot (Parameters::Rn, Parameters::a0, Parameters::mp);
    const Schroddy sh(potTot, H);
    ElementEigenfunctions elEigf = nuclei.orbitalEigenfunction(sh, orbitals);
    NuclearDensity NDens;
    NDens.theoreticalDensity(elEigf, nuclei.getLevelDegeneration());
    NDens.sogDensity(qrParam, H);

    //std::ofstream fOut(outputPath + "refDensity.txt");
    fOut.open(outputPath + "refInitialDensity.txt");
    fOut << NDens.getTheoreticalDensity();
    fOut.close();


   /* for(int i = 0; i < elEigf.size(); ++i)
    {
    	fOut.open(outputPath + "level" + i + "refInitialEigenFunction.txt");
        for (const auto& it : elEigf[i])
            fOut << it.first << "\t" << it.second << std::endl;
    }*/

    /*ElementEigenValues initialEigenvalues = nuclei.getLevelEigenvalue();
    fOut.open(outputPath + "refInitialEigenvalues.txt");
    fOut << initialEigenvalues;
    fOut.close();*/
    fOut.open(outputPath + "refSogDensity.txt");
    fOut << NDens.getSOGDensity();
    fOut.close();

    //Test Kohn-Sham inversion for initial harmonic potential
    fOut.open(outputPath + "refFirstKSPotential.txt");
    KohnShamInverse ksi(potTot, H);
    KohnShamInverse tempKsi = ksi;
    ksi.KSinverse(NDens, tempKsi);
    fOut << ksi.getKSPot();
    fOut.close();

    unsigned long loops = 0;
    while (!NDens.hasConverged()) //simulation loop
    //while (!ksi.hasConverged(tempKsi)) //simulation loop
    {
        const PotOut po(ksi, Parameters::mn, 0);
        const Schroddy sh_(po, H);
        elEigf = nuclei.orbitalEigenfunction(sh_, orbitals);
        NDens.theoreticalDensity(elEigf, nuclei.getLevelDegeneration());
        KohnShamInverse tempKsi_ = ksi;
        ksi.KSinverse(NDens, tempKsi_);

        ++loops;
        //if (loops%10 == 0)
            std::cerr << "Convergence distance: " << NDens.distanceToConvergence()
			//std::cerr << "Convergence distance: " << ksi.distanceToConvergence()
            << " after " << loops << " iterations." << std::endl;
    }

    fOut.open(outputPath + "refFinalDensity.txt");
    fOut << NDens.getTheoreticalDensity();;
    fOut.close();
    /*ElementEigenValues finalEigenvalues = nuclei.getLevelEigenvalue();
    fOut.open(outputPath + "refFinalEigenvalues.txt");
    fOut << finalEigenvalues;
    fOut.close();*/
    KSPotential finalPotential = ksi.getKSPot();
    fOut.open(outputPath + "refFinalPotential.txt");
    fOut << finalPotential;
    fOut.close();

    clock_t end = clock(); // Finish time
    double hours = (((double)(end - start))/CLOCKS_PER_SEC)/3600;

    std::cout << "CONVERGENCE IS DONE in: " << hours << " hours!" << " GREAT JOB!" << std::endl;
    std::cout << "Program executed successfully." << std::endl;
    return 0;
}
