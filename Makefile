CC = g++
CFLAGS = -O2 -std=c++11

GWSAT: GWSAT.cpp
	$(CC) $(CFLAGS) $< -o $@
