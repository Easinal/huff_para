CC = g++
CFLAGS = -O3 -mcx16 -march=native -std=c++17 -pthread -I/home/csmajs/zwan018/.parlaylib/usr/local/include -Wall -Wextra

all:	huff_para

huff_para:	huff_para.cpp
	$(CC) $(CFLAGS) huff_para.cpp -o huff_para
	
clean:
	rm -f huff_para
