FC=mpifort
SRC=./src
FLAGS=-march=native -frecursive -O2
LIBS=-llapack -lblas
TARGETS=$(SRC)/matsolve_serial.for

main:
	$(FC) -o main.exe $(SRC)/main.for $(TARGETS) $(FLAGS) $(LIBS)