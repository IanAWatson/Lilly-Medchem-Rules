# Location for the installation

.PHONY: supporting_libraries
.PHONY: Molecule
.PHONY: test
.PHONY: all

all: supporting_libraries Molecule

supporting_libraries:
	cd supporting_libraries && make install

Molecule:
	cd Molecule && make install

install: supporting_libraries Molecule
	@echo "No separate installation required" >&2

test: supporting_libraries Molecule
	cd test && ./dotest.sh

clean:
	cd supporting_libraries && make clean
	cd Molecule && make clean

uninstall:
	cd supporting_libraries && make uninstall
	cd Molecule && make uninstall
