# Makefile

CXX = g++
CXXFLAGS = -std=c++11 -O2 -Wall -Wextra -Wno-unused-parameter -I/usr/local/include
LDFLAGS = -L/usr/local/lib -lgmpxx -lgmp -lntl -lpthread -lcryptopp 
TARGET = hydra
SRC = hydra.cpp

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC) $(LDFLAGS)

clean:
	rm -f $(TARGET)
