FFLAGS = -O
OBJS  =  nleq2.o wnorm.o linalg_nleq2.o zibconst.o zibsec.o zibmon.o

nleq2: main_nleq2.o $(OBJS)
	$(FC) $(DFLAGS) main_nleq2.o $(OBJS) -o nleq2

nleq2-easy: main_nleq2_easy.o $(OBJS) nleq2e.o
	$(FC) $(DFLAGS) main_nleq2_easy.o $(OBJS) nleq2e.o -o nleq2-easy

clean:
	rm -f main_nleq2.o main_nleq2_easy.o nleq2e.o $(OBJS) nleq2 nleq2-easy
