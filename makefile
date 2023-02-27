GCC = g++-11
CPPFLAGS = -std=c++17 -O3
LINKFLAGS = $(CPPFLAGS)

report: clean_data latex clean_tex clean_cpp

prepare_data: build
	./main 1.0 1.0 1e-7
	./main 1.0 1.0 1e-9
	./main 1.0 1.0 1e-11 print
	./main 1.0 5.0 1e-7
	./main 1.0 5.0 1e-9
	./main 1.0 5.0 1e-11 print
	./main 0.01 1.0 1e-7
	./main 0.01 1.0 1e-9
	./main 0.01 1.0 1e-11 print
	echo "0.01,5.0,-,-,-,-" >>stats.csv

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