INCLUDE_DIRS = 
LIB_DIRS = 
CC = gcc
CPP = g++

CDEFS=
CFLAGS= -g -Wall -fopenmp $(INCLUDE_DIRS) $(CDEFS)
CFLAGS2= -O0 -g -Wall -fopenmp $(INCLUDE_DIRS) $(CDEFS)
LIBS=

PRODUCT= lorenz_openmp

HFILES=
CFILES= lorenz_openmp.c
CPPFILES=

SRCS= ${HFILES} ${CFILES}
OBJS= ${CFILES:.c=.o}

all:	${PRODUCT}

clean:
	-rm -f *.o *.NEW *~
	-rm -f ${PRODUCT} ${DERIVED} ${GARBAGE}

lorenz_openmp:	lorenz_openmp.c
	$(CC) $(CFLAGS) lorenz_openmp.c -o lorenz -lm
