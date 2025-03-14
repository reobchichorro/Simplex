FLAGS = -std=c++17 -g -O3

all : main.out #clear.o

clear.o:
	rm *.o

# main.out: main.o Utils.o Defs.o Alloc.o SubAlloc.o Situation.o GreedyLS.o MH.o
# 	g++ main.o Terrain.o Utils.o GuardType.o AllocGPos.o SubAlloc.o Situation.o GreedyLS.o ILS.o GA.o -o main.out $(FLAGS)

# MH.o: GreedyLS.o Situation.o Alloc.o Defs.o Utils.o
# 	g++ -c MH/ILS.cpp $(FLAGS)
# 	g++ -c MH/GA.cpp $(FLAGS)

# GreedyLS.o: Situation.o Alloc.o Defs.o Utils.o
# 	g++ -c MH/GreedyLS.cpp $(FLAGS)

# Situation.o: SubAlloc.o Alloc.o Defs.o Utils.o
# 	g++ -c Domain/Situation.cpp $(FLAGS)

# SubAlloc.o: Alloc.o Defs.o Utils.o
# 	g++ -c Domain/SubAlloc.cpp $(FLAGS)

# Alloc.o: Defs.o Utils.o
# 	g++ -c Domain/AllocGPos.cpp $(FLAGS)

# Defs.o: Utils.o
# 	g++ -c Domain/Terrain.cpp $(FLAGS)
# 	g++ -c Domain/GuardType.cpp $(FLAGS)

Utils.o:
	g++ -c Utils.cpp $(FLAGS)

main.o: main.cpp
	g++ -c main.cpp $(FLAGS)

clean:
	rm -r -f *.o *.out
