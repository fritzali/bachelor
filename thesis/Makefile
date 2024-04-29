all: main

main: FORCE | build
	latexmk \
		--lualatex \
		--output-directory=build \
		--interaction=batchmode \
		--halt-on-error \
	main.tex

build:
	mkdir -p build

clean:
	rm -rf build

FORCE:

.PHONY: all clean
