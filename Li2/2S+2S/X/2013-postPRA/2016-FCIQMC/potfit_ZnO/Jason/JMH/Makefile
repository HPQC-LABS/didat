#FFLAGS = -g -C -u -fullwarn
#FFLAGS = -O2 -C -mips4 -64 -r12000
#FFLAGS = -O2 -mips4 -64 -r12000
FFLAGS = -g -O2
OBJECTS = dspotfit.o masses.o readata.o readpotw.o vgen-betay.o alf.o cdjoel.o schrq.o nllssrr.o qrod.o fununc.o nmchng.o dyidpj-scatln-betay.o dvddp.o scatln.o bes.o wridata.o fixfreq.o dvcalc.o dampF.o prepot.o

fit: $(OBJECTS)
	f77 $(FFLAGS) $(OBJECTS) -L/central/lib -lnag -o dpot.x

$(OBJECTS): arrsizes.h

clean:
	rm $(OBJECTS)
