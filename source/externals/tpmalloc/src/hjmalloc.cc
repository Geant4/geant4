#include <iostream>
#include <stdio.h>
#include <time.h>
#include "tls.hh"
#include "tpmalloc/mymalloc.h"


G4ThreadLocal int mallocswitch = 0;
G4ThreadLocal mspace tlms = 0;

G4ThreadLocal long int alloccount = 0;
G4ThreadLocal long int freecount = 0;

G4ThreadLocal double alloctime = 0.0;
G4ThreadLocal double freetime = 0.0;

#define HJMALLOCDEBUG 0
#define HJMALLOCPROFILE 0

//#define TLMSSIZE 512
#define TLMSSIZE 256
//#define TLMSSIZE 128
//#define TLMSSIZE 32
//#define TLMSSIZE 16

void turnontpmalloc() {mallocswitch = 1;}
void turnofftpmalloc() {mallocswitch = 0;}
int  tpmallocflag() {return mallocswitch;}

void *malloc (size_t __size)
{
  #if HJMALLOCPROFILE
  alloccount++;
  #endif

  #if HJMALLOCDEBUG
  std::cout << "in malloc of mymalloc" << std::endl;
  #endif

  if (mallocswitch == 0)
    return mymalloc(__size);
  else
  {
    #if HJMALLOCDEBUG
    std::cout << "mspace in malloc of mymalloc" << std::endl;
    #endif
    if (tlms == 0) tlms = create_mspace(TLMSSIZE * 1024 * 1024, 0);

    #if HJMALLOCPROFILE
    double timed;
    struct timespec ts1, ts2;
    clock_gettime(CLOCK_REALTIME, &ts1);
    #endif

    void *tempmalloc = mspace_malloc(tlms, __size); 

    #if HJMALLOCPROFILE
    clock_gettime(CLOCK_REALTIME, &ts2);
    timed = (ts2.tv_sec - ts1.tv_sec) * 1000000000 + (ts2.tv_nsec - ts1.tv_nsec);
    alloctime = alloctime + timed;
    #endif

    return tempmalloc;
  }
}

void *calloc (size_t __nmemb, size_t __size)
{
  #if HJMALLOCPROFILE
  alloccount++;
  #endif

  #if HJMALLOCDEBUG
  std::cout << "in calloc of mymalloc" << std::endl;
  #endif
  if (mallocswitch == 0)
    return mycalloc(__nmemb, __size);
  else
  {
    #if HJMALLOCDEBUG
    std::cout << "mspace in alloc of mymalloc" << std::endl;
    #endif
    if (tlms == 0) tlms = create_mspace(TLMSSIZE * 1024 * 1024, 0);

    #if HJMALLOCPROFILE
    double timed;
    struct timespec ts1, ts2;
    clock_gettime(CLOCK_REALTIME, &ts1);
    #endif

    void *tempmalloc = mspace_calloc(tlms, __nmemb, __size);

    #if HJMALLOCPROFILE
    clock_gettime(CLOCK_REALTIME, &ts2);
    timed = (ts2.tv_sec - ts1.tv_sec) * 1000000000 + (ts2.tv_nsec - ts1.tv_nsec);
    alloctime = alloctime + timed;
    #endif

    return tempmalloc;
  }
}

void *realloc (void *__ptr, size_t __size)
{
  #if HJMALLOCPROFILE
  alloccount++;
  freecount++;
  #endif

  #if HJMALLOCDEBUG
  std::cout << "in realloc of mymalloc" << std::endl;
  #endif
  if (mallocswitch == 0)
    return myrealloc(__ptr, __size);
  else
  {
    #if HJMALLOCDEBUG
    std::cout << "mspace in realloc of mymalloc" << std::endl;
    #endif
    if (tlms == 0) tlms = create_mspace(TLMSSIZE * 1024 * 1024, 0);

    #if HJMALLOCPROFILE
    double timed;
    struct timespec ts1, ts2;
    clock_gettime(CLOCK_REALTIME, &ts1);
    #endif

    void *tempmalloc = mspace_realloc(tlms, __ptr, __size);

    #if HJMALLOCPROFILE
    clock_gettime(CLOCK_REALTIME, &ts2);
    timed = (ts2.tv_sec - ts1.tv_sec) * 1000000000 + (ts2.tv_nsec - ts1.tv_nsec);
    alloctime = alloctime + timed;
    #endif

    return tempmalloc;
  }
}

void free (void *__ptr)
{
  #if HJMALLOCPROFILE
  freecount++;
  #endif

  #if HJMALLOCDEBUG
  std::cout << "in free of mymalloc" << std::endl;
  #endif
  if (mallocswitch == 0)
    myfree(__ptr);
  else
  {
    #if HJMALLOCDEBUG
    std::cout << "mspace in free of mymalloc" << std::endl;
    #endif

    #if HJMALLOCPROFILE
    double timed;
    struct timespec ts1, ts2;
    clock_gettime(CLOCK_REALTIME, &ts1);
    #endif

    mspace_free(tlms, __ptr);

    #if HJMALLOCPROFILE
    clock_gettime(CLOCK_REALTIME, &ts2);
    timed = (ts2.tv_sec - ts1.tv_sec) * 1000000000 + (ts2.tv_nsec - ts1.tv_nsec);
    freetime = freetime + timed;
    #endif

  }
}

void cfree (void *__ptr)
{
  #if HJMALLOCPROFILE
  freecount++;
  #endif

  #if HJMALLOCDEBUG
  std::cout << "in cfree of mymalloc" << std::endl;
  #endif
  if (mallocswitch == 0)
    cfree(__ptr);
  else
  {
    #if HJMALLOCDEBUG
    std::cout << "mspace in cfree of mymalloc" << std::endl;
    #endif

    #if HJMALLOCPROFILE
    double timed;
    struct timespec ts1, ts2;
    clock_gettime(CLOCK_REALTIME, &ts1);
    #endif

    mspace_free(tlms, __ptr);

    #if HJMALLOCPROFILE
    clock_gettime(CLOCK_REALTIME, &ts2);
    timed = (ts2.tv_sec - ts1.tv_sec) * 1000000000 + (ts2.tv_nsec - ts1.tv_nsec);
    freetime = freetime + timed;
    #endif

  }
}

void *valloc (size_t __size)
{
  #if HJMALLOCPROFILE
  alloccount++;
  #endif

  #if HJMALLOCDEBUG
  std::cout << "in valloc of mymalloc" << std::endl;
  #endif
  if (mallocswitch == 0)
    return myvalloc(__size);
  else
  {
    #if HJMALLOCDEBUG
    std::cout << "mspace in valloc of mymalloc" << std::endl;
    #endif
    if (tlms == 0) tlms = create_mspace(TLMSSIZE * 1024 * 1024, 0);

    #if HJMALLOCPROFILE
    double timed;
    struct timespec ts1, ts2;
    clock_gettime(CLOCK_REALTIME, &ts1);
    #endif

    void *tempmalloc = mspace_malloc(tlms, __size);

    #if HJMALLOCPROFILE
    clock_gettime(CLOCK_REALTIME, &ts2);
    timed = (ts2.tv_sec - ts1.tv_sec) * 1000000000 + (ts2.tv_nsec - ts1.tv_nsec);
    alloctime = alloctime + timed;
    #endif

    return tempmalloc;
  }
}
