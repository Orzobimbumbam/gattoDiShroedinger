#include "theirMain.hpp"

void theirMain()
{
    const double H = 0.1;
    //mkdir("Outputs", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH); // Create outputs folder
    std::string inputPath = "../Inputs/";
    std::string outputPath = "Outputs/";
    clock_t start = clock(); // Start time

    //Load the matrix of quantum numbers for each state from "orbitals.txt"
    std::vector <std::vector <unsigned int>> orbitals (36);
    std::ifstream in (inputPath + "orbitals.txt");
    readMatrix(orbitals, in, false);
    in.close();

    in.open(inputPath + "40Ca.txt");
    std::vector<std::vector<double>> qrParam(12);
    readMatrix(qrParam, in, false);
    in.close();

    Element nuclei(orbitals);

    //Test theoretical and sog density for initial harmonic potential
    //const TotPot potTot (Parameters::mn);
    //const WSaxPot potTot (Parameters::Rn, Parameters::a0, Parameters::mp);
    const HOPot potTot(Parameters::mn); //default is ground state
    //const TestPot potTot(Parameters::mn);

    std::map<double, double> inPotential;
    const unsigned long NSteps = std::abs(Parameters::x_fin - Parameters::x_in)/H;
    double rad = Parameters::x_in;
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
    NDens.theoreticalDensity(elEigf, nuclei.getLevelDegeneration());
    NDens.benchmarkDensity(qrParam, H);

    fOut.open(outputPath + "refInitialDensity.txt");
    fOut << NDens.getTheoreticalDensity();
    fOut.close();
    fOut.open(outputPath + "refSogDensity.txt");
    fOut << NDens.getBenchmarkDensity();
    fOut.close();
    ElementEigenValues initialEigenvalues = nuclei.getLevelEigenvalue();
    fOut.open(outputPath + "refInitialEigenvalues.txt");
    writeMatrix(initialEigenvalues, fOut, false);
    fOut.close();
    fOut.open(outputPath + "refInitialEigenfunctions.txt");
    writeElementEigenfunctions(elEigf, fOut);
    fOut.close();
    /*
    NuclearDensityWithMC NDens_;
    NDens_.theoreticalDensity(elEigf, nuclei.getLevelDegeneration());
    std::vector<std::vector<double>> mcMatrix(NDens_.getTheoreticalDensity().size());
    std::ifstream fIn(outputPath + "refInitialDensity.txt");
    readMatrix(mcMatrix, fIn, false);
    NDens_.benchmarkDensity(mcMatrix);
     */
    //Test Kohn-Sham inversion for initial harmonic potential
    fOut.open(outputPath + "refFirstKSPotential.txt");
    KohnShamInverse ksi(potTot, H);
    KohnShamInverse tempKsi = ksi;
    //ksi.KSinverseWithLB(NDens, tempKsi);
    ksi.KSinverseWithJW(NDens, tempKsi);
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
        //ksi.KSinverseWithLB(NDens, tempKsi_);
        ksi.KSinverseWithJW(NDens, tempKsi);

        ++loops;
        //if (loops%10 == 0)
            std::cerr << "Convergence distance: " << NDens.distanceToConvergence() << " with epsilon "
            << NDens.epsilon() << " after " << loops << " iterations. " << std::endl;
    }

    fOut.open(outputPath + "refFinalDensity.txt");
    fOut << NDens.getTheoreticalDensity();;
    fOut.close();
    ElementEigenValues finalEigenvalues = nuclei.getLevelEigenvalue();
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
    
    return;
}