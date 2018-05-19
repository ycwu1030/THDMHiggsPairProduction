
LT = /Users/ycwu/Library/Mathematica/Applications/LoopTools/x86_64-Darwin/lib
THDMC = /Volumes/working-HDD/Users/teddy/workingspace/EWPM/2HDMC-1.7.0.bac/lib
LHAPDF = /Volumes/working-HDD/Users/teddy/workingspace/Useful-Package/LHAPDF/6.1.6/lib

SRCDIR := src
INCDIR := include
OBJDIR := obj
CXX = $(LT)/../bin/f++
FFLAG = -I$(LT)/../include -I$(INCDIR) -I$(THDMC)/../src `gsl-config --cflags` -I$(LHAPDF)/../include
FLIBS = -L$(LT) -looptools `gsl-config --libs` -L$(THDMC) -l2HDMC -lHB -lHS -L/usr/local/gfortran/lib -lgfortran -L$(LHAPDF) -lLHAPDF
RPATH = -Xlinker -rpath $(LHAPDF)

INC = $(wildcard $(INCDIR)/*.h)
SRC = $(wildcard $(SRCDIR)/*.cpp)
OBJ = $(patsubst $(SRCDIR)/%.cpp, $(OBJDIR)/%.o, $(SRC))

.PRECIOUS: $(OBJ)
.PHONY: all clean MKOBJDIR

all: HiggsPairProd.x MKOBJDIR

%.x: %.cpp $(OBJ)
	$(CXX) $< $(FFLAG) $(OBJ) $(FLIBS) $(RPATH) -o $@

$(OBJDIR)/pcubature.o: $(SRCDIR)/pcubature.cpp $(INCDIR)/cubature.h
	$(CXX) $(FFLAG) -c $< -o $@

$(OBJDIR)/hcubature.o: $(SRCDIR)/hcubature.cpp $(INCDIR)/cubature.h
	$(CXX) $(FFLAG) -c $< -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp $(INC)
	$(CXX) $(FFLAG) -c $< -o $@

MKOBJDIR:
	@mkdir -p $(OBJDIR)

clean:
	rm -f *.x
	rm -f $(OBJDIR)/*.o

