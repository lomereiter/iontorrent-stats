DMD=ldmd2
D_FLAGS=-O3 -release -inline -version=serial
D_INCLUDES=-Ilib/BioD -Isrc

all:
	mkdir -p build
	rdmd --build-only --compiler=$(DMD) $(D_FLAGS) $(D_INCLUDES) -ofbuild/iontorrent-stats src/D/collectstats.d 

.PHONY: clean

clean:
	rm -rf build/
