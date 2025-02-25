FLAGS=-O3 -mfpmath=sse -fstack-protector-all -g -W -Wall -Wextra -Wunused -Wcast-align -Werror -pedantic -pedantic-errors  -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-suggest-attribute=format -c

all: main.o reader.o matrix.o solve.o
	g++ main.o reader.o matrix.o solve.o

main.o: main.cpp func.h
	g++ $(FLAGS) main.cpp

reader.o: reader.cpp func.h
	g++ $(FLAGS) reader.cpp

matrix.o: matrix.cpp func.h
	g++ $(FLAGS) matrix.cpp

solve.o: solve.cpp func.h
	g++ $(FLAGS) solve.cpp

# results.o: results.cpp func.h
# 	g++ $(FLAGS) results.cpp

# basic_functions.o:basic_functions.cpp func.h
# 	g++ $(FLAGS) basic_functions.cpp

clean:
	rm -f *.out *.o *.gch
