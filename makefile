GCC = g++-10
CPPFLAGS = -std=c++17 -O3
LINKFLAGS = $(CPPFLAGS)

report: clean_data latex clean_tex clean_cpp

a1.0_b1.0: build
	./main 1.0 1.0

a1.0_b5.0: build
	./main 1.0 5.0

a0.01_b1.0: build
	./main 0.01 1.0

a0.01_b5.0: build
	./main 0.01 5.0

build: main.o
	${GCC} $? -o main $(LINKFLAGS)

%.o: %.cpp
	${GCC} -c $< -o $@ $(CPPFLAGS)

latex: a1.0_b1.0 a1.0_b5.0 a0.01_b1.0
	pdflatex -interaction=nonstopmode report > /dev/null
	pdflatex -interaction=nonstopmode report > /dev/null

clean_tex:
	rm -f *.aux *.log

clean_cpp:
	rm -f *.o main

clean_data:
	rm -f *.csv

clean: clean_tex clean_cpp clean_data
	rm -f report.pdf