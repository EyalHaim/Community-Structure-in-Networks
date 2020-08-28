FLAGS = -ansi -Wall -Wextra -Werror -pedantic-errors
LIBS = -lm

all: cluster
clean: 
		rm -rf *.o cluster

cluster: main.o spmat.o clustering.o divisionUtils.o IO.o list.o modularityMaximization.o 
		gcc main.o spmat.o clustering.o divisionUtils.o IO.o list.o modularityMaximization.o -o cluster $(LIBS)

main.o: main.c IO.h spmat.h clustering.h
		gcc $(FLAGS) -c main.c

spmat.o: spmat.c spmat.h errorHandler.h divisionUtils.h
		gcc $(FLAGS) -c spmat.c

clustering.o: clustering.c list.h spmat.h clustering.h errorHandler.h divisionUtils.h modularityMaximization.h
		gcc $(FLAGS) -c clustering.c

divisionUtils.o: divisionUtils.c divisionUtils.h errorHandler.h
		gcc $(FLAGS) -c divisionUtils.c

IO.o: IO.c spmat.h errorHandler.h IO.h
		gcc $(FLAGS) -c IO.c

list.o: list.c list.h
		gcc $(FLAGS) -c list.c
 
modularityMaximization.o: modularityMaximization.c modularityMaximization.h divisionUtils.h  errorHandler.h
		gcc $(FLAGS) -c modularityMaximization.c
  	



