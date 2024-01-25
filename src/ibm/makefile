EXE=patch_size.exe
CPP=main.cpp patch_size.cpp individual.cpp patch.cpp
HPP=patch_size.hpp individual.hpp patch.hpp
CXX=g++
CXXFLAGS=-Wall -ggdb -std=c++20

$(EXE) : $(CPP) $(HPP)
	$(CXX) $(CXXFLAGS) -o $(EXE) $(CPP)

clean: 
	rm -rf $(EXE)
