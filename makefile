GSL_LIBS=$(shell gsl-config --libs) 
CFLAGS=$(gsl-config --cflags) -Wall -O3
INCDIRS=-I/usr/local/include
LIBS=-L/usr/local/lib $GSL_LIBS
hka: hka.c hka.h
	gcc -o hka hka.c $(GSL_LIBS) $(INCDIRS) $(CFLAGS)


