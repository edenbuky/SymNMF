# -*- MakeFile -*-

mysymnmf.cpython-36m-x86_64-linux-gnu.so: symnmfmodule.c symnmf.o setup.py
	python3 setup.py build_ext --inplace


symnmf.o: symnmf.c symnmf.h
	gcc -o symnmf symnmf.c -lm
