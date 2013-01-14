struct MemorySegment {
  void *start;
  void *end;
};

static struct MemorySegment memorySegments[4096];
static char libname[4096][1024];
static int numOfSegments = 0;
