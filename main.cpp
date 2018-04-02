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

    in.open(inputPath + "48Ca.txt");
    std::vector<std::vector<double>> qrParam(12);
    readMatrix(qrParam, in, false);

    Element nuclei(orbitals);

    //Test theoretical and sog density for initial harmonic potential
    //const HOPot pot(Parameters::mn); //default is ground state
    const WSaxPot pot(Parameters::mn); //default is ground state
    const Schroddy sh(pot, H);
    ElementEigenfunctions elEigf = nuclei.orbitalEigenfunction(sh, orbitals);
    NuclearDensity NDens;
    NDens.theoreticalDensity(elEigf, nuclei.getLevelDegeneration());
    NDens.sogDensity(qrParam, H);

    std::ofstream fOut(outputPath + "refDensity.txt");
    fOut << NDens.getTheoreticalDensity();

    fOut.close();
    fOut.open(outputPath + "refSogDensity.txt");
    fOut << NDens.getSOGDensity();

    //Test Kohn-Sham inversion for initial harmonic potential
    fOut.close();
    fOut.open(outputPath + "refFirstKSPotential.txt");
    KohnShamInverse ksi(pot, H);
    KohnShamInverse tempKsi = ksi;
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
        KohnShamInverse tempKsi_ = ksi;
        ksi.KSinverse(NDens, tempKsi_);

        ++loops;
        //if (loops%10 == 0)
            std::cerr << "Convergence distance: " << NDens.distanceToConvergence()
            << " after " << loops << " iterations." << std::endl;
    }


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
