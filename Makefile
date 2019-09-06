FORT = gfortran
CFLAGS = -c
LIBS = -llapack -lblas

all: delaunayLP data
	./generate
	./delaunayLP deldata.txt

delaunayLP: delaunayLPtest.f90 dualsimplex.o
	$(FORT) delaunayLPtest.f90 dualsimplex.o $(LIBS) -o delaunayLP

dualsimplex.o: dualsimplex.f90
	$(FORT) $(CFLAGS) dualsimplex.f90 -o dualsimplex.o

data: generate_data.f90
	$(FORT) generate_data.f90 -o generate

clean:
	rm -f *.o *.mod delaunayLP generate
