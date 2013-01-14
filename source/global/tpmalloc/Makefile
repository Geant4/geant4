.PHONY: all clean

all: message libtpmallocstub.so libmymalloc.so

libtpmallocstub.so: tpmallocstub.c
	@gcc -fPIC -O3 -c tpmallocstub.c
	@gcc -shared -O3 -o libtpmallocstub.so tpmallocstub.o

libmymalloc.so: hjmalloc.c mymalloc.c
	@gcc -fPIC -O3 -c -o hjmalloc.o hjmalloc.c
	@gcc -fPIC -O3 -c -o mymalloc.o mymalloc.c
	@gcc -shared -O3 -o libmymalloc.so mymalloc.o hjmalloc.o

message:
	@echo "Building tpmalloc libraries ..."
clean:
	@rm -rf *.o *.so
