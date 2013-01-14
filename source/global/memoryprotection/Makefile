all: libProtectedMemory.so tracer

PROTECTMEMORY_OBJS = ProtectedMemory.o

libProtectedMemory.so: $(PROTECTMEMORY_OBJS)
	gcc -shared -o $@ $(PROTECTMEMORY_OBJS)

ProtectedMemory.o: ProtectedMemory.c
	gcc -fPIC -c $<

TRACER_OBJS = tracer.o filter.o

tracer: $(TRACER_OBJS)
	gcc -g -O0 -o $@ -lpthread $(TRACER_OBJS)

tracer.o: tracer.c
	gcc -g -c $<

filter.o: filter.c
	gcc -g -c $<

clean:
	rm -f *.o *.d *.so *~ sharedupdate tracer memoryhotpot.out.* gdbcommands.*

