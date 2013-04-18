#ifndef ProtectedMemory_h
#define ProtectedMemory_h 1

struct MemorySegment {
  void *start;
  void *end;
};

static struct MemorySegment memorySegments[4096];
static int numOfSegments = 0;

#endif
