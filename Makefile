#HCOMPILER := $(HOME)/Apps/gcc-14/bin/g++
HCOMPILER := g++

FLAGS := -lquadmath -static-libstdc++ -static-libgcc

SRC:=./src
BIN:=./bin

$(shell mkdir -p $(BIN))

BINARIES=$(BIN)/example_test
SRCS = example_test.cpp PowerSeries.cpp

all: $(BINARIES)

$(BIN)/example_test: $(BIN)/example_test.o $(BIN)/PowerSeries.o
	$(HCOMPILER) -o $@ $^ $(FLAGS)

$(BIN)/%.o : $(SRC)/%.cpp
	$(HCOMPILER) -c -o $@ $^ $(FLAGS)

clean:
	@$(RM) -rf $(BIN)

cleanall:
	@$(RM) -rf $(BIN)
