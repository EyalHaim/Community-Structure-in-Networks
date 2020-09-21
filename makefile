FLAGS = -ansi -Wall -Wextra -Werror -pedantic-errors
LIBS = -lm

all: cluster
clean: 
		rm -rf *.o cluster

cluster: main.o graph.o division.o IO.o list.o  
		gcc main.o graph.o division.o IO.o list.o -o cluster $(LIBS)

main.o: main.c IO.h graph.h division.h
		gcc $(FLAGS) -c main.c

graph.o: graph.c graph.h errorHandler.h
		gcc $(FLAGS) -c graph.c

division.o: division.c list.h graph.h division.h errorHandler.h 
		gcc $(FLAGS) -c division.c
		
IO.o: IO.c graph.h errorHandler.h IO.h
		gcc $(FLAGS) -c IO.c

list.o: list.c list.h
		gcc $(FLAGS) -c list.c
 
  	



