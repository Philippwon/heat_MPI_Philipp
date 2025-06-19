# Intel compiler
CC =  icc
CFLAGS = -Ofast -xCORE-AVX512 -fno-alias
MPICC = mpicc




all: heat

heat : heat.o input.o misc.o timing.o relax_gauss.o relax_jacobi.o
	$(CC) $(CFLAGS) -o $@ $+ -lm 

%.o : %.c heat.h timing.h input.h
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -f *.o heat *~ *.ppm

remake : clean all






all_MPI: heat_MPI

heat_MPI : heat_MPI.o input.o misc.o timing.o relax_gauss.o relax_jacobi.o
	$(MPICC) $(CFLAGS) -o $@ $+ -lm 

%.o : %.c heat.h timing.h input.h
	$(MPICC) $(CFLAGS) -c -o $@ $<

clean_MPI:
	rm -f *.o heat_MPI *~ *.ppm

remake_MPI : clean all