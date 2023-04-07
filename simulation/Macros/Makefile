# MakeFile

CPPFLAGS = \
	-I$(includedir) \
	-I$(OFFLINE_MAIN)/include \
	-I$(ROOTSYS)/include \
	-I$(shell `root-config --incdir`)

LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64 
  

BINS = $(patsubst %.C, %, $(wildcard *.C))


default : $(BINS)

$(BINS): % : %.C 
	@echo -n "Building $@ ... "
	$(CXX) $(CCFLAGS) $< -o $@ -w $(LDFLAGS) $(shell fastjet-config --libs --cxxflags) $(shell root-config --cflags --glibs) $(shell pythia8-config --cxxflags --libs) -lpythia8 -ldl
	@echo "Done"


clean: 
	rm -f *.o	$(BINS)

jobclean:
	rm *.o* *.e*
