This is an example of how to write a Python Extention in C.

The FASTAOP_core.c contains the C function definitions, and the FASTAOP.c contains the
Pyhton interface.

All the methods defined in FASTAOP_methods in the FASTAOP.c file are available from
Python.

Compilation:
------------

	# python setup.py build


Installation:
------------

	# sudo python setup.py install


Running tests:
--------------

	# cp build/LIBRARY_DIR/LIB.so .
	# python test.py

See test.py for more info on how to use the extension.
