CC := gcc
CFLAGS := -DPOLYBENCH_TIME -O3 -march=native
LOADLIBES := -lrt

SRC := atax.c lu.c jacobi-2d.c
OBJ := $(SRC:.c=.o) polybench.o
BIN := $(SRC:.c=)

.PHONY: all clean

all : $(BIN)

atax : polybench.o atax.o
	$(LINK.c) $^ $(LDFLAGS) $(LOADLIBES) -o $@
jacobi-2d : polybench.o jacobi-2d.o
	$(LINK.c) $^ $(LDFLAGS) $(LOADLIBES) -o $@
lu : polybench.o lu.o
	$(LINK.c) $^ $(LDFLAGS) $(LOADLIBES) -o $@

clean :
	-rm -vf $(BIN) $(OBJ)
