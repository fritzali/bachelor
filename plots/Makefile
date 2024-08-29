all: main present

build:
	mkdir -p build

clean:
	rm -rf build

FORCE:

.PHONY: all clean

main: code Makefile main.py style.rc header.tex | build
	TEXINPUTS=$$(pwd): MATPLOTLIBRC=style.rc python main.py
	touch build/main

present: code Makefile present.py option.rc design.tex | build
	TEXINPUTS=$$(pwd): MATPLOTLIBRC=option.rc python present.py
	touch build/present
