FILES = mergeArcBPM


all: $(FILES)


$(FILES): % : %.cc
	@echo "Building $@ ... "

	g++ `root-config --libs --cflags` -o $@ $@.cc
#	g++ -o $@ $@.cc

	@echo "Done"

