
LT = /Users/ycwu/Library/Mathematica/Applications/LoopTools/x86_64-Darwin/lib
THDMC = /Volumes/working-HDD/Users/teddy/workingspace/EWPM/2HDMC-1.7.0.bac/lib
#LT = /Users/mac/work/LoopTools/x86_64-Darwin/lib

SRCDIR := src
INCDIR := include
OBJDIR := obj
CXX = $(LT)/../bin/f++
FFLAG = -I$(LT)/../include -I$(INCDIR) -I$(THDMC)/../src `gsl-config --cflags`
FLIBS = -L$(LT) -looptools `gsl-config --libs` -L$(THDMC) -l2HDMC -lHB -lHS -L/usr/local/gfortran/lib -lgfortran

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

