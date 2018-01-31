//ROOTFLAGS = `root-config --cflags`
//LIBROOT = `root-config --libs`

empidensity.o: empidensity.cpp
	g++ -c empidensity.cpp $(ROOTFLAGS) $(LIBROOT)

initpot.o: initpot.cpp
	g++ -c initpot.cpp $(ROOTFLAGS) $(LIBROOT)

kohn-sham.o: kohn-sham.cpp
	g++ -c kohn-sham.cpp $(ROOTFLAGS) $(LIBROOT)

schroddy.o: schroddy.cpp
	g++ -c schroddy.cpp $(ROOTFLAGS) $(LIBROOT)

main.o: main.cpp
	g++ -c main.cpp $(ROOTFLAGS) $(LIBROOT)

Build: empidensity.o initpot.o kohn-sham.o schroddy.o main.o
	g++ empidensity.o initpot.o kohn-sham.o schroddy.o main.o -o exec.x $(ROOTFLAGS) $(LIBROOT)

Execute: Execute exec.x
	./exec.x

