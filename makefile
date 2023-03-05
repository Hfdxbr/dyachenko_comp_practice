GCC = g++-11
CPPFLAGS = -std=c++17 -O3
LINKFLAGS = $(CPPFLAGS)

report: clean_data latex clean_tex clean_cpp

prepare_data: build
	./main 50 50 1.0
	./main 100 100 1.0
	./main 200 200 1.0
	./main 50 50 0.0
	./main 100 100 0.0
	./main 200 200 0.0

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