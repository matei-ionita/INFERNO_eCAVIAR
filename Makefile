CC=g++
DIC=$(PWD)
LDFLAGS= -O1 -I $(DIC)/armadillo/include/ -DARMA_DONT_USE_WRAPPER -llapack -lblas -lgsl

SRC=main.cpp Model.cpp utils.cpp
OBJ=$(SRC:.cpp = .o)

INFERNO_eCAVIAR: $(SRC)
	$(CC) $(SRC) -o INFERNO_eCAVIAR $(LDFLAGS)

clean:
	rm INFERNO_eCAVIAR
