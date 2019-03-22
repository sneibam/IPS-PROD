all: src

src:
	$(MAKE) -C src
test:
	$(MAKE) -C src test
doc :
	doxygen
format:
	astyle --options=astyle.conf src/*.cpp,*.h
.PHONY: all src test clean doc
clean:
	$(MAKE) -C src clean
	
