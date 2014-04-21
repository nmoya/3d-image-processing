LIB=./lib
INCLUDE=./include
SRC=./src
OBJ=./obj
BIN=./bin

#FLAGS= -g  -O0 -Wall -D _DEBUG -Wno-unused-result -fPIC -std=gnu99 -pedantic
FLAGS=-O3 -Wall -fPIC -std=gnu99 -pedantic -I /System/Library/Frameworks/vecLib.framework/Headers/

libmo815-3dvis: $(LIB)/libmo815-3dvis.a
	echo "libmo815-3dvis.a built..."

$(LIB)/libmo815-3dvis.a: \
$(OBJ)/common.o \
$(OBJ)/adjacency.o \
$(OBJ)/kernel.o \
$(OBJ)/image.o \
$(OBJ)/matrix.o 
	ar csr $(LIB)/libmo815-3dvis.a \
$(OBJ)/common.o \
$(OBJ)/adjacency.o \
$(OBJ)/kernel.o \
$(OBJ)/image.o \
$(OBJ)/matrix.o \

$(OBJ)/common.o: $(SRC)/common.c
	gcc $(FLAGS) -c $(SRC)/common.c -I$(INCLUDE) -o $(OBJ)/common.o 

$(OBJ)/adjacency.o: $(SRC)/adjacency.c
	gcc $(FLAGS) -c $(SRC)/adjacency.c -I$(INCLUDE) -o $(OBJ)/adjacency.o 

$(OBJ)/kernel.o: $(SRC)/kernel.c
	gcc $(FLAGS) -c $(SRC)/kernel.c -I$(INCLUDE) -o $(OBJ)/kernel.o 

$(OBJ)/image.o: $(SRC)/image.c
	gcc $(FLAGS) -c $(SRC)/image.c -I$(INCLUDE) -o $(OBJ)/image.o 

$(OBJ)/matrix.o: $(SRC)/matrix.c
	gcc $(FLAGS) -c $(SRC)/matrix.c -I$(INCLUDE) -o $(OBJ)/matrix.o 

clean: 
	rm $(LIB)/lib*.a; rm $(OBJ)/*.o; rm $(BIN)/* 	





