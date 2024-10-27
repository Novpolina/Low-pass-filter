CFLAGS = -Wall -Wextra

all: furie.o
	@g++ $(CFLAGS) furie.o -o furie 

signal.o: firie.cpp
	@g++ -c furie.cpp 

run_signal: furie.o
	./furie
	python3 signal_lpf.py

run_spectrum: furie.o
	./furie
	python3 spectrum.py

clean:
	del *.o *.exe *.txt
