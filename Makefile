all : DNA check

check : test_DNA

DNA : 
	python3 setup.py build
	cp build/lib*/*.so .

test_DNA : 
	python3 -m pytest

clean :
	@rm -rf __pycache__ build .pytest_cache *.so

.PHONY : clean all