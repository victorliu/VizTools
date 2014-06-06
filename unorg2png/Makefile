#CFLAGS = -O0 -ggdb
CFLAGS = -O3

unorg2png: deltri.o lodepng.o read_data_1d.o line_reader.o predicates.o main.cpp
	$(CXX) $(CFLAGS) -I.. deltri.o lodepng.o read_data_1d.o line_reader.o predicates.o main.cpp -o unorg2png
deltri.o: deltri.c
	$(CC) -c $(CFLAGS) deltri.c -o deltri.o
lodepng.o: lodepng.cpp
	$(CXX) -c $(CFLAGS) lodepng.cpp -o lodepng.o
read_data_1d.o: read_data_1d.c
	$(CC) -c $(CFLAGS) -I. read_data_1d.c -o read_data_1d.o
line_reader.o: line_reader.c
	$(CC) -c $(CFLAGS) -I. line_reader.c -o line_reader.o
predicates.o: predicates.c
	$(CC) -c $(CFLAGS) predicates.c -o predicates.o

clean:
	rm -f *.o unorg2png
