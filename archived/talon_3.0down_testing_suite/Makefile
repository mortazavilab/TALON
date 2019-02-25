test :
	if [ -d "scratch" ]; then rm -Rf scratch; fi
	mkdir scratch;
	pytest;

unit :
	pytest -v -m unit;

integration :
	if [ -d "scratch" ]; then rm -Rf scratch; fi
	mkdir scratch;
	pytest -v -m integration;

clean :
	if [ -d "scratch" ]; then rm -Rf scratch; fi
	mkdir scratch;
