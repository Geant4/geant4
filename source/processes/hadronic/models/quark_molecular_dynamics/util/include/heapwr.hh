//
// Define macros for heap debug
//

#ifndef __WRAPPER__

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#ifdef CHECK_HEAP
#define HEAPWRAP
#else
#undef HEAPWRAP
#endif

#ifdef IS_SGI
#undef HEAPWRAP
#endif

#ifdef HEAPWRAP
#define __HDB__
#undef __NHDB__
#undef __CLEAR__
#undef NEW
#endif
#ifndef HEAPWRAP
#undef __HDB__
#undef __NHDB__
#define __CLEAR__
#endif
// #undef __WRAPPER__

#ifdef __HDB__
   #define NEW new(__FILE__, __LINE__)
   #define DELETE(x) ((delete (x)),(x) = 0)
#endif
      
#ifdef __NHDB__
   #define NEW new
   #define DELETE(x) ((delete (x)), (x) = 0)
#endif
            
#ifdef __CLEAR__
   #define NEW new
   #define DELETE(x) ((delete (x)))
#endif

#define NAMELENGTH 20
//
// Basisklasse Wrapper f"ur die "Uberpr"ufung des Heaps
//
struct Wrapper {
    // Prolog
    Wrapper *pNext;
    int length;
    char fileName[NAMELENGTH];
    unsigned int lineNumber;
    char prologue;
    int :0;      // User-Daten sollen auf eine word boundary fallen
    // User data
    char data[1];
    // Epilog folgt hier
};


// Prototyp-Deklarationen

void werror(Wrapper*, char*, int = 0);

//
// "Uberladung der new-Operatoren
//

#ifndef __CLEAR__
void* operator new(size_t);
void* operator new(size_t, char*, int);
#ifdef IS_GCC
void* operator new[](size_t);
void* operator new[](size_t, char*, int);
#endif
#endif
void* allocateFn(size_t sizeBlock, char *pFile, int lineNo);

//
// "Uberladung der delete-Operatoren
//
#ifndef __CLEAR__
void operator delete(void*);
#endif
void deleteFn(void*);
void displayAllocated();

#define __WRAPPER__
#endif

