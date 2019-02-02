test :
	set -e;
	make clean;
	python build_test_databases.py;
	pytest;

unit :
	set -e;
	make clean;
	pytest -v -m unit;

integration :
	set -e;
	make clean;
	pytest -v -m integration;

clean :
	set -e; 	
	if [ -d "scratch" ]; then rm -Rf scratch; fi
	mkdir scratch;
