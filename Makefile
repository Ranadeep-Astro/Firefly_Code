
UNAME = $(shell uname)
ifeq ($(UNAME),Linux)
H55 = /usr/lib/x86_64-linux-gnu/hdf5/serial
endif
ifeq ($(UNAME),Darwin)
H55 = /opt/local
endif

CC = mpicc
FLAGS = -O3 -Wall -g

INC = -I$(H55)/include
LIB = -L$(H55)/lib

OBJ = main.o misc.o rad.o h5io.o readpar.o

default: glow

%.o: %.c header.h structs.h
	$(CC) $(FLAGS) $(INC) -c $<

 
glow: $(OBJ) header.h structs.h
	$(CC) $(FLAGS) $(LIB) -o glow $(OBJ) -lhdf5 -lm

clean:
	rm -f *.o glow
