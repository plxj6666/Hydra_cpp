# Makefile

CXX = g++
# 添加 -g 选项生成调试信息，并将 -O2 优化降低到 -O0，以便更好地进行调试
CXXFLAGS = -std=c++11 -g -O0 -Wall -Wextra -Wno-unused-parameter  -I../include -I/usr/local/include \
			-I$(XKCP)/lib/high/Keccak \
            -I$(XKCP)/lib/common \
            -I$(XKCP)/lib/low/KeccakHash \
		    -I$(XKCP)/bin/generic64/libXKCP.so.headers
LDFLAGS = -L/usr/local/lib -L$(XKCP)/bin/generic64 \
		 -lgmpxx -lgmp -lntl -lpthread -lcryptopp -lXKCP \
		  -Wl,-rpath=$(XKCP)/bin/generic64
XKCP = ../../test_module/XKCP
TARGET = hydra_change
SRC = hydra_change.cpp

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC) $(LDFLAGS)

clean:
	rm -f $(TARGET)
