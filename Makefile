CC = g++
CFLAGS = -std=c++17

TARGET = simulace

all: $(TARGET)

$(TARGET): $(TARGET).cpp
	$(CC) $(CFLAGS) -o $(TARGET) $(TARGET).cpp 

clean:
	$(RM) $(TARGET) *.plt *.png
