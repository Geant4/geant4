PROG   = X11Torus
CFLAGS = -Wall -O3
CC     = gcc
LIB    = -L/usr/lib/ -lm -L/usr/X11R6/lib/ -lm -lX11
OBJ    = xdemo.o main.o torus.o

$(PROG):	$(OBJ)
		$(CC) $(CFLAGS) -s $(LIB) $(OBJ) -o $(PROG)

xdemo.o:	xdemo.c xdemo.h
main.o:		main.c main.h xdemo.h

clean:
	rm -f *.o
