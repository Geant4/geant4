// see license file for original license.

#ifndef tools_glutess_memalloc
#define tools_glutess_memalloc

////////////////////////////////////////////////////////
/// inlined C code : ///////////////////////////////////
////////////////////////////////////////////////////////

#include <cstdlib>
#include <cstring>
namespace tools {
inline void* mem_alloc(size_t a_size) {return ::memset( ::malloc( a_size ), 0xa5, a_size);}
inline void* mem_realloc(void* a_ptr,size_t a_size){return ::realloc(a_ptr,a_size);}
inline void mem_free(void* a_ptr){::free(a_ptr);}
inline int mem_init( size_t /*maxFast*/) {return 1;}
}

#define memAlloc	tools::mem_alloc
#define memRealloc	tools::mem_realloc
#define memFree		tools::mem_free
#define memInit		tools::mem_init

#endif
