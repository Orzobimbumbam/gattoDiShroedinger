//ROOTFLAGS = `root-config --cflags`
//LIBROOT = `root-config --libs`

parameters.o: parameters.cpp
	g++ -c parameters.cpp $(ROOTFLAGS) $(LIBROOT)
	
IOUtils.o: IOUtils.cpp
	g++ -c IOUtils.cpp $(ROOTFLAGS) $(LIBROOT)

densities.o: densities.cpp
	g++ -c densities.cpp $(ROOTFLAGS) $(LIBROOT)
	
idensity.o: densities.cpp
	g++ -c idensity.cpp $(ROOTFLAGS) $(LIBROOT)

mcDensity.o: densities.cpp
	g++ -c mcDensity.cpp $(ROOTFLAGS) $(LIBROOT)

eigenfunction.o: eigenfunction.cpp
	g++ -c eigenfunction.cpp $(ROOTFLAGS) $(LIBROOT)

eigenvalues.o: eigenvalues.cpp
	g++ -c eigenvalues.cpp $(ROOTFLAGS) $(LIBROOT)

element.o: element.cpp
	g++ -c element.cpp $(ROOTFLAGS) $(LIBROOT)

initpot.o: initpot.cpp
	g++ -c initpot.cpp $(ROOTFLAGS) $(LIBROOT)

kohn-sham.o: kohn-sham.cpp
	g++ -c kohn-sham.cpp $(ROOTFLAGS) $(LIBROOT)

schroddy.o: schroddy.cpp
	g++ -c schroddy.cpp $(ROOTFLAGS) $(LIBROOT)

main.o: main.cpp
	g++ -c main.cpp $(ROOTFLAGS) $(LIBROOT)

Build: parameters.o IOUtils.o densities.o idensity.o mcDensity.o eigenfunction.o eigenvalues.o element.o initpot.o kohn-sham.o schroddy.o main.o
	g++ parameters.o IOUtils.o densities.o idensity.o mcDensity.o eigenfunction.o eigenvalues.o element.o initpot.o kohn-sham.o schroddy.o main.o -o exec.x $(ROOTFLAGS) $(LIBROOT)

Execute: Execute exec.x
	./exec.x

