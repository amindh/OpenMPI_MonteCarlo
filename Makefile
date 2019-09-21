OPTI=-O
#OPTI=-g -O0
OPT=O_c++

CXX=/usr/local/bin/mpicxx $(OPTI)
CC=/usr/local/bin/mpicc $(OPTI)

MPI_DIR=/usr/local
MPI_INC=-I$(MPI_DIR)/include
MPI_LIB=-L$(MPI_DIR)/lib
MPI_LIBS= -lmpi

LD=/usr/local/bin/mpicxx $(OPTI)

EXE= monte_carlo_para

C++SOURCES= monte_carlo_para

SOURCESFILES= $(C++SOURCES:=.cc)
OBJECTSFILES= $(C++SOURCES:=.o)

default : all

all : $(OBJECTSFILES)

	@ echo "Edition de lien"
	$(LD) -o $(EXE) ./*.o $(MPI_INC) $(MPI_LIB) $(MPI_LIBS) -ldl

%.o: %.c
	echo "Compilation de " $< " : "
	$(CXX) $(MPI_INC) -c $<




