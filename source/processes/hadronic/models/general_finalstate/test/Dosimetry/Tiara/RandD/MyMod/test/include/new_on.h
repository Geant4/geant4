#ifdef __DEBUG__

#ifdef __INIT_MEMORY__

#define new new(__FILE__,__LINE__)
#define delete (SetOwner(__FILE__,__LINE__),false)? SetOwner("",0) : delete

#define malloc(x) AllocateMemory(x,__FILE__,__LINE__,MEM_MALLOC)
#define calloc(x,y) AllocateMemory((x)*(y),__FILE__,__LINE__,MEM_CALLOC)
#define realloc(x,y) ReallocateMemory(x,y,__FILE__,__LINE__,MEM_REALLOC)
#define free(x) deAllocateMemory(x,__FILE__,__LINE__,MEM_FREE)

#endif

#endif
