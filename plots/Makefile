all: build/main

build:
	mkdir -p build

clean:
	rm -rf build

FORCE:

.PHONY: all clean

build/main: figure-config Makefile main.py header-matplotlib.tex | build
	TEXINPUTS=$$(pwd): MATPLOTLIBRC=figure-config python main.py
	touch build/main
