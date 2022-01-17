CC = g++
#CFLAGS = -O3 -mcx16 -march=native -std=c++17 -pthread -I/home/csmajs/zwan018/.parlaylib/usr/local/include -Wall -Wextra
CFLAGS = -O3 -mcx16 -march=native -std=c++17 -DPARLAY_CILKPLUS -fcilkplus -DCILK -I/home/csmajs/zwan018/.parlaylib/usr/local/include -Wall -Wextra
#CFLAGS = -O3 -mcx16 -march=native -std=c++17 -DPARLAY_SEQUENTIAL -I/home/csmajs/zwan018/.parlaylib/usr/local/include -Wall -Wextra

all:	huff_para test

huff_para:	huff_para.cpp
	$(CC) $(CFLAGS) huff_para.cpp -o huff_para
	
test: test.cpp
	$(CC) $(CFLAGS) test.cpp -o test

clean:
	rm -f huff_para
