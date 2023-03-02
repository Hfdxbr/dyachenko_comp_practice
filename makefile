GCC = g++-11
CPPFLAGS = -std=c++17 -O3
LINKFLAGS = $(CPPFLAGS)

report: clean_data latex clean_tex clean_cpp

prepare_data: build
	./main 1e-7
	./main 1e-9
	./main 1e-11 print

build: main.o
	${GCC} $? -o main $(LINKFLAGS)

%.o: %.cpp
	${GCC} -c $< -o $@ $(CPPFLAGS)

latex: prepare_data
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