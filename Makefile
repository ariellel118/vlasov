CC = g++
LDFLAGS = `pkg-config --libs ibsimu-1.0.5new_solver`
CXXFLAGS = -Wall -g `pkg-config --cflags ibsimu-1.0.5new_solver`

vlasov2d: vlasov2d.o
	$(CC) -o vlasov2d vlasov2d.o $(LDFLAGS)

vlasov2d.o: vlasov2d.cpp
	$(CC) -c -o vlasov2d.o vlasov2d.cpp $(CXXFLAGS)

clean:	
	$(RM) *~ *.o vlasov2d