#ifdef __DEBUG__

#ifndef __INIT_MEMORY__
#define __INIT_MEMORY__
#endif

#ifdef __INIT_MEMORY__
#ifndef __MEMORY_INITIALIZED__
#define __MEMORY_INITIALIZED__

#include "new_off.h"
#include "stdlib.h"
void InitializeManager();
void* AllocateMemory(size_t size,char* FileName,unsigned LineNumber,char Method);
void* ReallocateMemory(void* pAddress,size_t size,char* FileName,unsigned LineNumber,char Method);
void deAllocateMemory(void* pvAddress,char* FileName,unsigned lineNumber,char Method);
void SetOwner(char* FileName,unsigned LineNumber);

class InitMemManager{
public:
  InitMemManager(){InitializeManager();}
  ~InitMemManager(){;};
};

#ifndef __MMAN_EXISTS__
#define __MMAN_EXISTS__
InitMemManager mMan;
#endif

#define MEM_UNKNOWN      0
#define MEM_NEW          1
#define MEM_NEW_ARRAY    2
#define MEM_DELETE       3
#define MEM_DELETE_ARRAY 4
#define MEM_MALLOC       5
#define MEM_CALLOC       6
#define MEM_REALLOC      7
#define MEM_FREE         8

inline void* operator new(size_t size,char* FileName,unsigned LineNumber){
  return AllocateMemory(size,FileName,LineNumber,MEM_NEW);
};

inline void* operator new[](size_t size,char* FileName,unsigned LineNumber){
  return AllocateMemory(size,FileName,LineNumber,MEM_NEW_ARRAY);
}

inline void operator delete(void* pAddress){
  if(pAddress != NULL)
    deAllocateMemory(pAddress,"",0,MEM_DELETE);
}
inline void operator delete[](void* pAddress){
  if(pAddress != NULL)
    deAllocateMemory(pAddress,"",0,MEM_DELETE_ARRAY);
}
#endif
#include "new_on.h"
#endif

#endif
