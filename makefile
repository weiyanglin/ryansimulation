LIBS = -L/usr/X11R6/lib -lX11
FLAGS = -O2 -Wno-write-strings
acoustics1D: acoustics1D.o
	g++ *.o $(LIBS) $(FLAGS) -o acoustics1D
acoustics1D.o: acoustics1D.cpp
	g++ acoustics1D.cpp $(LIBS) $(FLAGS) -c
clean:
	rm acoustics1D *.o
