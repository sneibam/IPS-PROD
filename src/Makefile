CC = g++
CFLAGS = -std=c++11 -Wall
TESTGEN = cxxtestgen
TESTFLAGS = --error-printer
TESTS = UnitTest.h
TESTSCPP = Poly.cpp Basis.cpp Miscellaneous.cpp Density.cpp
OBJSTEST = Poly.o Basis.o Miscellaneous.o Density.o
DEPS= Miscellaneous.h
TARGET = main
TARGETTEST = test
OBJS = Poly.cpp Miscellaneous.cpp Basis.cpp Density.cpp main.cpp

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $^ -o $@ -larmadillo -llapack -lblas
	./$(TARGET)

%.o: %.cpp %.h
	$(CC) $(CFLAGS) $< -c -o $@

$(TARGETTEST): $(TESTS)
	$(CC) $(CFLAGS) -c $(TESTSCPP) -larmadillo
	$(TESTGEN) $(TESTFLAGS) -o test.cpp $(TESTS)
	$(CC) $(CFLAGS) -o $@  test.cpp $(OBJSTEST) -larmadillo -llapack -lblas
	./$@

.PHONY: test clean

clean:
	rm -f $(TARGET) $(TARGETTEST) *.o t.cpp test.cpp *.orig