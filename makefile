GCC = g++-11
CPPFLAGS = -std=c++17 -O3
LINKFLAGS = $(CPPFLAGS)

report: clean_data latex clean_tex clean_cpp

prepare_data: build
	./main 51 51 0.1 0.1
	./main 101 101 0.1 0.1
	./main 201 201 0.1 0.1
	./main 51 51 1.0 0.1
	./main 101 101 1.0 0.1
	./main 201 201 1.0 0.1
	./main 51 51 1.0 1.0
	./main 101 101 1.0 1.0
	./main 201 201 1.0 1.0
	/bin/python3 relation.py

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
	rm -f report.pdf report.synctex.gz report.fls report.fdb*