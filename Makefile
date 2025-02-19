all: lorentz_taylor.exe
CXXFLAGS=-std=c++17 -O3 -DNDEBUG -march=native -ffast-math -I $(BOOST_ROOT)/include

ctaylor/cjacobian.h:
	git submodule init ctaylor
	git submodule update --init --recursive ctaylor

clean:
	rm -f lorenz_taylor.exe

lorentz_taylor.exe:lorenz_taylor.cpp ctaylor/cjacobian.h
	$(CXX) $(CXXFLAGS) -o $@ lorenz_taylor.cpp
