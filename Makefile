all: main

CXX = clang++
override CXXFLAGS += -g -Wall -Werror -std=c++20 -fsanitize=address,undefined

SRCS = integer.cpp
LIBS = *.h

main: $(SRCS) $(LIBS) main.cpp
	$(CXX) $(CXXFLAGS) $(SRCS) main.cpp -o "$@"

bdz2: $(SRCS) $(LIBS) bdz2.cpp
	$(CXX) $(CXXFLAGS) $(SRCS) bdz2.cpp -o "$@"

main-debug: $(SRCS) $(LIBS) main.cpp
	$(CXX) $(CXXFLAGS) -O0 $(SRCS) main.cpp -o "$@"

clean:
	rm -rf *.dSYM main main-debug bdz2
