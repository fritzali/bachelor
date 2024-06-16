all: build/main

build:
	mkdir -p build

clean:
	rm -rf build

FORCE:

.PHONY: all clean

build/main: code Makefile main.py style.rc header.tex | build
	TEXINPUTS=$$(pwd): MATPLOTLIBRC=style.rc python main.py
	touch build/main
