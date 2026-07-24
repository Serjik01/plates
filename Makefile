CXX      ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra
TARGET    = plate
HEADERS   = plate.hpp linalg.hpp quad.hpp basis.hpp navier.hpp galerkin.hpp kantorovich.hpp output.hpp

$(TARGET): main.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) main.cpp -o $(TARGET)

run: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(TARGET)
	rm -rf out

.PHONY: run clean
