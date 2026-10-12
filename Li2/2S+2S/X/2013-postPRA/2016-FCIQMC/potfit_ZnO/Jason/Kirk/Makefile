#FFLAGS = -O2
FFLAGS = -g -C -fullwarn
#FFLAGS = -g -C -u -fullwarn
OBJECTS = hecokp3d.o kpcoef.o

bcont: $(OBJECTS)
        f77 $(FFLAGS) $(OBJECTS) -o kirk.x
