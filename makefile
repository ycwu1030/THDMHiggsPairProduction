
LT = /Users/ycwu/Library/Mathematica/Applications/LoopTools/x86_64-Darwin/lib
#LT = /Users/mac/work/LoopTools/x86_64-Darwin/lib

SRCDIR := src
INCDIR := include
OBJDIR := obj
CXX = $(LT)/../bin/f++
FFLAG = -I$(LT)/../include -I$(INCDIR) `gsl-config --cflags`
FLIBS = -L$(LT) -looptools `gsl-config --libs`

SRC = $(wildcard $(SRCDIR)/*.cpp)
OBJ = $(patsubst $(SRCDIR)/%.cpp, $(OBJDIR)/%.o, $(SRC))

all: HiggsPairProd.x

%.x:%.cpp $(OBJ)
	$(CXX) $(FFLAG) -o $@ $< $(OBJ) $(FLIBS)

MKOBJDIR:
	@mkdir -p $(OBJDIR)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp MKOBJDIR
	$(CXX) $(FFLAG) -c $< -o $@

.PHONY: clean

clean:
	rm -f *.x
	rm -f $(OBJDIR)/*.o

