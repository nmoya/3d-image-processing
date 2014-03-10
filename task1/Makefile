#FLAGS= -g -O0 -Wall -fopenmp -pthread -std=gnu99 -pedantic -lm
FLAGS= -O3 -Wall -fopenmp -pthread -std=gnu99 -pedantic -lm

$@.c: 
	gcc $(FLAGS) nmImage.c $@.c -o $@;
