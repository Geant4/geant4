#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <libgen.h> /* for basename() */
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/wait.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <errno.h>
#include <assert.h>
//#include <ucontext.h>
#include <sys/mman.h>
#include <signal.h>
#include <execinfo.h>
#include <stdlib.h>
#include <sys/syscall.h>
#include <pthread.h>

#include "tls.hh"

static pid_t gettid()
{
  pid_t mytid = syscall(SYS_gettid);
  return mytid;
}


static char *startPosition = 0;
static char *currentPosition = 0;
static char *nextPosition = 0;
static char *endPosition = 0;

void* AllocateInProtectedHeap(size_t size)
{
  if (startPosition != 0)
  {
    currentPosition = nextPosition;
    nextPosition = &nextPosition[size];
    if (nextPosition >= endPosition)
    {
      nextPosition = currentPosition;
      return malloc(size);
    }
    return currentPosition;
  }
  else
  {
    return malloc(size);
  }
}

static void AllocateProtectedHeap(size_t totalSize, size_t pagesize)
{

  if (startPosition != 0) return;

  startPosition = (char *) malloc(pagesize * 2 + totalSize);
  currentPosition = startPosition;
  nextPosition = currentPosition + pagesize;
  endPosition = nextPosition + totalSize;
}


#include "memoryprotection/UDSignals.h"

static int readFile(char *filename, char *buf, int maxSize);
static char *readhex(char *string, long *value);
static char *readMapsLine(char *string, void **start, void **end, int *writable, char *filename);
static void *getPageNumber(void *theAddress);
static void setuphandler0();
static void setuphandler1();
static void setuphandler2();
// static void raisesignal();
static void SetUpAllProtect();
static void UnSetUpAllProtect();
// static void protectAll();
// static void unProtectAll();
// static void noProtectAll();


static size_t pagesize;
static int logpagesize;
// static int stdOutFile;
// static int theSignal;
// static void *theAddress;

#include "memoryprotection/MemorySegment.h"

static char application[1024];

//#define PROTECTEDMEMORYDEBUG 1

static int phase = 0;

static void *getPageNumber(void *address) {return (void *) (((long) address >> logpagesize) << logpagesize);}

static G4ThreadLocal  void *trace[100];
static int outputFile = -1;
static int trace_size;

void PhaseChange()
{
  int tmpoutputFile;
  tmpoutputFile = outputFile;

  phase++;
  pid_t myselfpid = getpid();
  char outputfilename[1024];  
  outputfilename[0] = 0;
  sprintf(outputfilename, "./memoryhotpot.out.%d.%d", phase, myselfpid);

  outputFile = open(outputfilename, O_CREAT|O_RDWR, S_IRWXU);
  close(tmpoutputFile);
}

static void record_stackframe() {
  trace_size = 0;
  trace_size = backtrace(trace, 100);
  backtrace_symbols_fd(trace, trace_size, outputFile);
}

static void show_stackframe() {
  trace_size = 0;
  trace_size = backtrace(trace, 100);
  backtrace_symbols_fd(trace, trace_size, STDOUT_FILENO);
}

static int isCalled = 0;

static void handler0(int, siginfo_t *, void *)
{
  #ifdef PROTECTEDMEMORYDEBUG
  printf("In handler 0 p1.\n");
  #endif
  if (isCalled == 0)
  {
    #ifdef PROTECTEDMEMORYDEBUG
    printf("In handler 0 p2.\n");
    #endif
    isCalled = 1;
    SetUpAllProtect();
    #ifdef PROTECTEDMEMORYDEBUG
    printf("In handler 0 p3.\n");
    #endif
    assert( 0 == raise(SIGUSR2));
  }
}

static void handler1(int, siginfo_t *, void *)
{
  #ifdef PROTECTEDMEMORYDEBUG
  printf("In handler 1 p1.\n");
  #endif

  UnSetUpAllProtect();

  //  show_stackframe();

  record_stackframe();

  #ifdef PROTECTEDMEMORYDEBUG
  printf("In handler 1 p2.\n");
  #endif

  assert( 0 == raise(SIGUSR2));
}

// mprotect_map and mprotect_addr should be the statically defined arrays
// this function protect only one page, it is not correct because
// we do not know how large the memory is in one assemble instruction
static void handler2(int, siginfo_t *, void *)
{
  #ifdef PROTECTEDMEMORYDEBUG
  printf("In handler 2 p1.\n");
  #endif

  SetUpAllProtect();
}

void BuildProtectedMemory(size_t totalHeap)
{
  pagesize = getpagesize();

  pid_t myselfpid = getpid();
  char outputfilename[1024];  
  outputfilename[0] = 0;
  sprintf(outputfilename, "./memoryhotpot.out.%d.%d", phase, myselfpid);
  outputFile = open(outputfilename, O_CREAT|O_RDWR, S_IRWXU);

  show_stackframe();

  #ifdef PROTECTEDMEMORYDEBUG
  printf("Page size: %d\n", pagesize);
  #endif

  long tmp;
  for (logpagesize = -1, tmp = pagesize; tmp > 0; logpagesize++)
    tmp = tmp >> 1;

  #ifdef PROTECTEDMEMORYDEBUG
  printf("log(page size): %d\n", logpagesize);
  #endif

  if (totalHeap > 0)
  {
    AllocateProtectedHeap(totalHeap, pagesize);

    memorySegments[numOfSegments].start = getPageNumber(nextPosition);
    memorySegments[numOfSegments].end = getPageNumber(endPosition); 
    numOfSegments++;
  }

  //force to load backtrace library
  //  show_stackframe();

  //open output file
  //  assert(-1 != (stdOutFile = open("./memoryhotspot.out", O_CREAT, S_IRWXU)));

  char buf1[1024];
  char buf2[51200];
  char *tmp2 = buf2;
  static char filename2[1024];
  void *start;
  void *end;
  int writable;

  if (-1 == readFile((char *) "/proc/self/cmdline", buf1, 1024)) {
    perror("readFile: /proc/self/cmdline");
    abort();
  }
  memset(application, 0, 1024);
  strcpy(application, buf1);

  //  cout << "command line: " << application << endl;

  if (-1 == readFile((char *) "/proc/self/maps", buf2, 51200)) {
    perror("readFile: /proc/self/maps");
    abort();
  }
  static char lastfilename[1024];
  do {
    tmp2 = readMapsLine(tmp2, &start, &end, &writable, filename2);
    if (strlen(filename2) > 0)
    {
      lastfilename[0] = 0;
      strcpy(lastfilename, filename2); 
    }
    else
    {
      strcpy(filename2, lastfilename); 
    }

    //    if (writable == 1 && (long int) start < 0x700000000000 && 
    //        strncmp(filename2, "/lib64/", 7) != 0 && strncmp(filename2, "/opt/", 5) != 0 &&
    //        ((long int) nextPosition < (long int) start || (long int) nextPosition > (long int ) end) && 
    //        ((long int) endPosition < (long int) start || (long int) endPosition > (long int ) end))
    if (writable ==  1)
    {
      //ignore 
      if (strstr(filename2, "ld") != 0)
      {
        std::cout << "Not protected: " << start << ", " << end << ", "
                                       << writable << ", " << filename2
                  << std::endl;
      }
      else if (strstr(filename2, "vdso") != 0)
      {
        std::cout << "Not protected: " << start << ", " << end << ", "
                                       << writable << ", " << filename2
                  << std::endl;
      }
      else if (strstr(filename2, "stack") != 0)
      {
        std::cout << "Not protected: " << start << ", " << end << ", "
                                       << writable << ", " << filename2
                  << std::endl;
      }
      else if (strstr(filename2, "libc") != 0)
      {
        std::cout << "Not protected: " << start << ", " << end << ", "
                                       << writable << ", " << filename2
                  << std::endl;
      }
      else if (strstr(filename2, "libProtectedMemory.so") != 0)
      {
        std::cout << "Not protected: " << start << ", " << end << ", "
                                       << writable << ", " << filename2
                  << std::endl;
      }
      else if (strstr(filename2, "libParRunManagerMP.so") != 0)
      {
        std::cout << "Not protected: " << start << ", " << end << ", "
                                       << writable << ", " << filename2
                  << std::endl;
      }
      else if (strstr(filename2, "libmymalloc.so") != 0)
      {
        std::cout << "Not protected: " << start << ", " << end << ", "
                                       << writable << ", " << filename2
                  << std::endl;
      }
      else if (strstr(filename2, "libtpmallocstub.so") != 0)
      {
        std::cout << "Not protected: " << start << ", " << end << ", "
                                       << writable << ", " << filename2
                  << std::endl;
      }
      else if (strstr(filename2, "heap") != 0)
      {
        static int Heapcount = 0;
        //Only one is the data segment follows with [heap] [thread private heap, empty, thread local storage]*          
        Heapcount++;
        if (Heapcount == 1)
        {
          std::cout << "Protected: " << numOfSegments << ", " << start << ", " << end << ", "
                                     << writable << ", " << filename2
                    << std::endl;
          memorySegments[numOfSegments].start = start;
          memorySegments[numOfSegments].end = end;
          numOfSegments++;
        }
        else
        {
          std::cout << "Not protected: " << start << ", " << end << ", "
                                       << writable << ", " << filename2
                    << std::endl;
        }
      }
      else if (strstr(filename2, "ParA01") != 0)
      {
        static int Parmaincount = 0;
	//Only one is the data segment follows with [heap] [thread private heap, empty, thread local storage]*
        Parmaincount++;
        printf("Determine: %d\n", Parmaincount);
	//        if (Parmaincount == 1 || Parmaincount == 2)
	{ 
          std::cout << "Protected: " << numOfSegments << ", " << start << ", " << end << ", "
                                     << writable << ", " << filename2
                    << std::endl;
          memorySegments[numOfSegments].start = start;
          memorySegments[numOfSegments].end = end;
          numOfSegments++;
	}
      }
      else if (strstr(filename2, "geant4") != 0)
      {
        std::cout << "Protected: " << numOfSegments << ", " << start << ", " << end << ", "
                                   << writable << ", " << filename2
                  << std::endl;
        memorySegments[numOfSegments].start = start;
        memorySegments[numOfSegments].end = end;
        numOfSegments++;
      }

    }
    else
    {
      std::cout << "Not protected: " << start << ", " << end << ", "
                                     << writable << ", " << filename2
                << std::endl;
    }

  } while ((*tmp2) != '\0' && numOfSegments < 4096);
}

static int readFile(char *filename, char *buf, int maxSize) {
  int num;
  int fd = open(filename, O_RDONLY);
  char *buf2 = buf;

  errno = 0;
  do {
    num = 0;
    do {
      buf2 += num;
      num = read(fd, buf2, 1000);
    } while (num != 0 && num != -1 && buf2 - buf < maxSize - 1000);
  } while (errno == EINTR || errno == EAGAIN);
  buf2 += num;
  *buf2 = '\0';
  *(buf2 + 1) = '\0'; 
  close(fd);
  return (errno == 0 ? 0 : -1);
}

static char *readhex (char *string, long *value) {
  char c;
  *value = 0;
  while (1) {
    c = *(string++);
    if ((c >= '0') && (c <= '9')) c -= '0';
    else if ((c >= 'a') && (c <= 'f')) c -= 'a' - 10;
    else if ((c >= 'A') && (c <= 'F')) c -= 'A' - 10;
    else break;
    *value = *value * 16 + c;
  }
  return (char *) --string;
}

static char *readMapsLine(char *string, void **start, void **end, int *writable, char *filename) {
  char *str;

  string = readhex(string, (long *) start) + 1;
  string = readhex(string, (long *) end);
  *writable = (string[2] == 'w');
  for (str = string; *str != '\n' && *str != '\0'; str++);
  *(str++) = '\0';
  for (str = string; (*str != '/' && *str != '[') && *str != '\0'; str++);

  filename[0] = 0;
  strcpy(filename, str);

  while (*str != '\0')
    str++;
  return str+1;
}

static void SetUpAllProtect()
{
  int i;
  for(i = 0 ; i < numOfSegments ; i++)
  {
    struct MemorySegment *tmpSegment = memorySegments + i;
    #ifdef PROTECTEDMEMORYDEBUG
    printf("Protect from 0x%lx to 0x%lx.\n", tmpSegment->start, tmpSegment->end);
    #endif

    if (-1 == mprotect(tmpSegment->start, (long) tmpSegment->end - (long) tmpSegment->start , PROT_READ))
      perror("mprotect error");

    #ifdef PROTECTEDMEMORYDEBUG
    printf("Finish to protect from 0x%lx to 0x%lx.\n", tmpSegment->start, tmpSegment->end);
    #endif
  }
}

static void UnSetUpAllProtect()
{
  int i;
  for(i = 0 ; i < numOfSegments ; i++)
  {
    struct MemorySegment *tmpSegment = memorySegments + i;

    #ifdef PROTECTEDMEMORYDEBUG
    printf("Remove protection from 0x%lx to 0x%lx.\n", tmpSegment->start, tmpSegment->end);
    #endif

    if (-1 == mprotect(tmpSegment->start, (long) tmpSegment->end - (long) tmpSegment->start , PROT_READ | PROT_WRITE))
      perror("mprotect error");
  }
}

#define MAXTHREADS 128
int numberOfThreads = 0;
pid_t threads[MAXTHREADS];
pid_t tracerpid;

//main thread
void addThread()
{
  numberOfThreads++;
  assert(numberOfThreads <= MAXTHREADS);
  threads[numberOfThreads-1] = gettid();
}


//forward declaration of a function used here
void setuphandlers();

//main thread
void starttracer()
{
  BuildProtectedMemory(0);
  setuphandlers();

  if ((tracerpid = fork()) == 0) //child process
  {
    char processid[10240];
    char tmptid[1024];
    processid[0] = 0;
    tmptid[0] = 0;
    sprintf(tmptid, "%d", numberOfThreads);
    strcat(processid, tmptid);
    strcat(processid, " ");
    int i;
    for (i = 0; i < numberOfThreads; i++)
    {
      tmptid[0] = 0;
      sprintf(tmptid, "%d", threads[i]);
      strcat(processid, tmptid);
      strcat(processid, " ");
    }
    execl("./memoryprotection/tracer", "./memoryprotection/tracer", application, processid, (char *)0);
    exit(0);
  }

}

//main thread
//There is only one worker thread allowed.
//This function is invoked after the worker thread exits.
void finishtracer()
{
  UnSetUpAllProtect();
  close(outputFile);
  assert(0 == raise(SIGUSR1));
  printf("Before wait for child process: %d\n", tracerpid);
  waitpid(tracerpid, NULL, 0);
  printf("After wait for child process: %d\n", tracerpid);
}

//worker thread
void StartDetection()
{
  printf("waiting for the tracer\n");
  while (isCalled == 0) usleep(100);
  printf("the tracer ready\n");
  assert(0 == raise(SIGUSR2));
}

//worker thread
void FinishDetection()
{
  //move to main thread
  //  UnSetUpAllProtect();

  assert(0 == raise(SIGUSR1));
  //For each thread, can not do this
}

/*
void starttracermain()
{
  //  pid_t currentID = getpid();
  pid_t currentID = gettid();
  if (fork()== 0) //child process
  {
    char processid[1024];
    processid[0] = 0;
    sprintf(processid, "%d", currentID);
    sleep(2);
    execl("./tracer", "./tracer", processid, application,  (char *)0);
    exit(0);
  }
  sleep(5);

}

void starttracerother()
{
  //  pid_t currentID = getpid();
  pid_t currentID = gettid();
  if (fork()== 0) //child process
  {
    char processid[1024];
    processid[0] = 0;
    sprintf(processid, "%d", currentID);
    sleep(2);
    execl("./tracer", "./tracer", processid, application, processid, (char *)0);
    exit(0);
  }
  sleep(5);

}
*/

void setuphandlers()
{
  setuphandler0();
  setuphandler1();
  setuphandler2();
}

static void setuphandler0()
{
  struct sigaction sa;
  sa.sa_sigaction = handler0;
  sigfillset(&sa.sa_mask);
  sa.sa_flags = SA_RESTART | SA_SIGINFO; /* Restart functions if
					    interrupted by handler */
  if (sigaction(SIGUSR2, &sa, NULL) == -1)
    perror("sigaction");  /* Handle error */;
}

static void setuphandler1()
{
  struct sigaction sa;
  sa.sa_sigaction = handler1;
  sigfillset(&sa.sa_mask);
  sa.sa_flags = SA_RESTART | SA_SIGINFO; /* Restart functions if
					    interrupted by handler */
  if (sigaction(SIGSEGV, &sa, NULL) == -1)
    perror("sigaction");  /* Handle error */;
}

static void setuphandler2()
{
  struct sigaction sa;

  sa.sa_sigaction = handler2;
  sigfillset(&sa.sa_mask);
  sa.sa_flags = SA_RESTART | SA_SIGINFO; /* Restart functions if
                                            interrupted by handler */
  if (sigaction(SIGUSR1, &sa, NULL) == -1)
    perror("sigaction");  /* Handle error */;
}

