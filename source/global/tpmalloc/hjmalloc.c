#include <stdio.h>
#include <time.h>
#include "mymalloc.h"


__thread int mallocswitch = 0;
__thread mspace tlms = 0;

__thread long int alloccount = 0;
__thread long int freecount = 0;

__thread double alloctime = 0.0;
__thread double freetime = 0.0;

#define HJMALLOCDEBUG 0
#define HJMALLOCPROFILE 0

//#define TLMSSIZE 512
#define TLMSSIZE 256
//#define TLMSSIZE 128
//#define TLMSSIZE 32
//#define TLMSSIZE 16

void turnontpmalloc() {mallocswitch = 1;};
void turnofftpmalloc() {mallocswitch = 0;};
int  tpmallocflag() {return mallocswitch;};

void *malloc (size_t __size)
{
  #if HJMALLOCPROFILE
  alloccount++;
  #endif

  #if HJMALLOCDEBUG
  printf("in malloc of mymalloc\n");
  #endif

  if (mallocswitch == 0)
    return mymalloc(__size);
  else
  {
    #if HJMALLOCDEBUG
    printf("mspace in malloc of mymalloc\n");
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
  printf("in calloc of mymalloc\n");
  #endif
  if (mallocswitch == 0)
    return mycalloc(__nmemb, __size);
  else
  {
    #if HJMALLOCDEBUG
    printf("mspace in alloc of mymalloc\n");
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
  printf("in realloc of mymalloc\n");
  #endif
  if (mallocswitch == 0)
    return myrealloc(__ptr, __size);
  else
  {
    #if HJMALLOCDEBUG
    printf("mspace in realloc of mymalloc\n");
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
  printf("in free of mymalloc\n");
  #endif
  if (mallocswitch == 0)
    myfree(__ptr);
  else
  {
    #if HJMALLOCDEBUG
    printf("mspace in free of mymalloc\n");
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
  printf("in cfree of mymalloc\n");
  #endif
  if (mallocswitch == 0)
    cfree(__ptr);
  else
  {
    #if HJMALLOCDEBUG
    printf("mspace in cfree of mymalloc\n");
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
  printf("in valloc of mymalloc\n");
  #endif
  if (mallocswitch == 0)
    return myvalloc(__size);
  else
  {
    #if HJMALLOCDEBUG
    printf("mspace in valloc of mymalloc\n");
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
