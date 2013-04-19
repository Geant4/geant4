#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <iostream>
#include <sys/types.h>
#include <sys/syscall.h>
#include <sys/ptrace.h>
#include <sys/wait.h>
#include <unistd.h>
#include <errno.h>
#include <sys/user.h>
#include <assert.h>
#include <sys/select.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>
#include <pthread.h>

#include "tls.hh"
#include "memoryprotection/filter.h"
#include "memoryprotection/UDSignals.h"

int fildes[2];
int status;
char cmd[1024];
int protect = 0;

#define MAXTHREADS 128

pid_t filter;

int numberOfThreads;
pid_t threads[MAXTHREADS];
pthread_t tids[MAXTHREADS];
static char libname[4096][1024];

#if defined(__MACH__) 
struct user_regs_struct
{
 long int ebx;
 long int ecx;
 long int edx;
 long int esi;
 long int edi;
 long int ebp;
 long int eax;
 long int xds;
 long int xes;
 long int xfs;
 long int xgs;
 long int orig_eax;
 long int eip;
 long int xcs;
 long int eflags;
 long int esp;
 long int xss;
 long int rip;
};
#define __WCLONE     0x80000000
#endif

G4ThreadLocal struct user_regs_struct regs;
G4ThreadLocal struct user_regs_struct regssinglestep;
G4ThreadLocal int detectionFlag = 0;

pthread_mutex_t segFaultHandlingLock = PTHREAD_MUTEX_INITIALIZER;

#include "memoryprotection/MemorySegment.h"

void FetchMap(pid_t pid)
{
  char fname[128], line[4096];
  sprintf(fname, "/proc/%d/maps", pid);
  FILE *fp = fopen(fname, "r");

  char *n1, *n2, *p;

  char lastlibname[1024];
    
  while (fgets(line, sizeof(line), fp) != NULL) {
    p = NULL;
    p = strstr(line, "/");
    
    if (p != NULL)
    {
      strcpy(lastlibname, p);
    }

    p = strstr(line, "[");

    if (p != NULL)
    {
      strcpy(lastlibname, p);
    }

    strcpy(libname[numOfSegments], lastlibname);

    n1 = strtok(line, "- ");
    n2 = strtok(NULL, "- ");
    memorySegments[numOfSegments].start = (void*)strtoul(n1, NULL, 16);
    memorySegments[numOfSegments].end = (void*)strtoul(n2, NULL, 16);

    numOfSegments++;
  }
  fclose(fp);

  int i;
  for (i = 0 ; i < numOfSegments ; i++)
  {
    std::cout << "Maps: " << memorySegments[i].start << ", "
              << memorySegments[i].end << ", " << libname[i] << std::endl;
  }
}

char *getLibName(long int address)
{
  int i;
  for (i = 0 ; i < numOfSegments ; i++)
  {
    if ((void*)address >= memorySegments[i].start && (void*)address <= memorySegments[i].end)
      return libname[i];
  }
  return NULL;
}

void processParameters(int, char *argv[])
{
  strcpy(cmd, argv[1]);
  printf("tracer: traced process executable file: %s:: %s %s\n", cmd, argv[1], argv[2]);

  char tmpnum[1024];
  char *current;
  current = argv[2];
  int i = 0;
  while (*current != ' ')
  {
    tmpnum[i] = *current;
    current++; i++;
  }

  tmpnum[i] = 0;
  current++;
  numberOfThreads = atoi(tmpnum);
  printf("tracer number of threads: %d.\n", numberOfThreads);

  int j;
  for (j = 0 ; j < numberOfThreads; j++)
  {
    i = 0;
    while (*current != ' ')
    {
      tmpnum[i] = *current;
      current++; i++;
    }
    tmpnum[i] = 0;
    current++;
    threads[j] = atoi(tmpnum);
    printf("tracer threads[%d] = %d\n", j, threads[j]);
  }
}

void test(pid_t pid, int wait_val)
{
  printf("tracer for %d WIFEXITED: %d\n", pid, WIFEXITED(wait_val));
  printf("tracer for %d WEXITSTATUS: %d\n", pid, WEXITSTATUS(wait_val));
  printf("tracer for %d WIFSIGNALED: %d\n", pid, WIFSIGNALED(wait_val));
  printf("tracer for %d WTERMSIG: %d\n", pid, WTERMSIG(wait_val));
  printf("tracer for %d WIFSTOPPED: %d\n", pid, WIFSTOPPED(wait_val));
  printf("tracer for %d WSTOPSIG: %d\n", pid, WSTOPSIG(wait_val));
  printf("tracer for %d WIFCONTINUED: %d\n", pid, WIFCONTINUED(wait_val));
}

static G4ThreadLocal int isMain = 1;

static int my_waitpid(pid_t pid, int *state, int flags)
{
  int ret;
  do
    {
      ret = waitpid(pid, state, flags);
    } while (ret == -1 && errno == EINTR);
  return ret;
}

static int my_wait(pid_t pid, int *wait_val)
{
  if (isMain)
    return my_waitpid(pid, wait_val, 0);
  else
    return my_waitpid(pid, wait_val, __WCLONE);
}

int SingleSteps(pid_t pid)
{
  int wait_val;           /*  traced's return value        */

#if defined(__linux__) 
  if (ptrace(PT_GETREGS, pid, 0, &regssinglestep) != 0)
  {
    printf("tracer for %d ptrace getregs", pid);
    exit(-1);
  }
#endif

#ifdef __x86_64__ 
  write(fildes[1], &(regssinglestep.rip), sizeof(regssinglestep.rip));
#else
  write(fildes[1], &(regssinglestep.eip), sizeof(regssinglestep.eip));
#endif

  if (ptrace(PT_STEP, pid, 0, 0) != 0)
  {
    printf("tracer for %d ptrace singlestep", pid);
    exit(-1);
  }
  assert(pid == my_wait(pid, &wait_val));
  if (wait_val == 1407) return 0;
  if (WSTOPSIG(wait_val) == SIGTERM)
  {
    printf("tracer for %d terminated signal\n", pid);
    exit(-1);
  }
  else if (wait_val == 0)
  {
    printf("tracer for %d finished signal\n", pid);
    exit(-1);
  }
  else
  {
    test(pid, wait_val);
    exit(-1);
  }
  return 0;
}

int traceloop(pid_t pid)
{
  int wait_val;

  if (ptrace(PT_ATTACH, pid, 0, 0) != 0)
  {
    printf("tracer for %d ptrace attach", pid);
    exit(-1);
  }
  assert(pid == my_wait(pid, &wait_val));

  printf("tracer for %d after attach\n", pid);

  if (isMain)
  {
    if (ptrace(PT_CONTINUE, pid, 0, SIGUSR2) != 0)
    {
      printf("tracer for %d ptrace continue with SIGUSR2", pid);
      exit(-1);
    }

    assert(pid == my_wait(pid, &wait_val));
    assert(SIGUSR2 == WSTOPSIG(wait_val));

    protect = 1;
    printf("tracer for %d after memory protection\n", pid);

    if (ptrace(PT_CONTINUE, pid, 0, 0) != 0)
    {
      printf("tracer for %d ptrace continue without SIGUSR2", pid);
      exit(-1);
    }

  }
  else
  {
    printf("tracer for %d Waiting for the memory is protected\n", pid);
    while (!protect) usleep(10);

    if (ptrace(PT_CONTINUE, pid, 0, 0) != 0)
    {
      printf("tracer for %d ptrace continue without SIGUSR2", pid);
      exit(-1);
    }

    printf("tracer for %d after continue\n", pid);

  }

  int tmp = 0;

  do
  {
    // printf("tracer for %d tracerloop\n", pid);
    // pid_t frompid = my_wait(pid, &wait_val);

    if (WSTOPSIG(wait_val) == SIGSEGV)
    {
      pthread_mutex_lock(&segFaultHandlingLock);

#if defined(__linux__) 
      if (ptrace(PT_GETREGS, pid, 0, &regssinglestep) != 0)
      {
        printf("tracer for %d ptrace getregs", pid);
      }
#endif

      if (detectionFlag && !isMain)
      {
#ifdef __x86_64__ 
        write(fildes[1], &(regssinglestep.rip), sizeof(regssinglestep.rip));
        printf("tracer for %d Address: %lx in %s", pid, regssinglestep.rip, getLibName(regssinglestep.rip));
#else
        write(fildes[1], &(regssinglestep.eip), sizeof(regssinglestep.eip));
        printf("tracer for %d Address: %lx in %s", pid, regssinglestep.eip, getLibName(regssinglestep.eip));
#endif
        char convertcmd[1024];
        memset(convertcmd, 0, 1024);
#ifdef __x86_64__ 
        sprintf(convertcmd, "addr2line -e %s 0x%lx", cmd, regssinglestep.rip);  
#else
        sprintf(convertcmd, "addr2line -e %s 0x%lx", cmd, regssinglestep.eip);  
#endif
        system(convertcmd);
      }

      if (ptrace(PT_CONTINUE, pid, 0, SIGSEGV) != 0)
      {
        printf("tracer for %d ptrace continue with SIGSEGV", pid);
	tmp = -6;
      }

      assert(pid == my_wait(pid, &wait_val));
      assert(SIGUSR2 == WSTOPSIG(wait_val));


      if (SIGUSR2 != WSTOPSIG(wait_val))
      {
        test(pid, wait_val);
        tmp = -3; 
      }
      else
      {
        SingleSteps(pid);
      //      counter++;
      //      printf("Counter: %d\n", counter);

        if (ptrace(PT_CONTINUE, pid, 0, SIGUSR1) != 0)
        {
          printf("tracer for %d ptrace continue with SIGUSR1\n", pid);
          tmp = -5;
        }
      }

      pthread_mutex_unlock(&segFaultHandlingLock);

    }
    else
    {
      if ( wait_val == 0 )
      {
	//        test(pid, wait_val);
        tmp = -1;
      }
      else if (WSTOPSIG(wait_val) == SIGUSR1)
      {
	//        test(pid, wait_val);
        detectionFlag = 0;
        if (ptrace(PT_DETACH, pid, 0, 0) != 0)
        {
          printf("tracer for %d ptrace detach\n", pid);
          tmp = -6;
        }
      }
      else if (WSTOPSIG(wait_val) == SIGUSR2)
      {
        detectionFlag = 1;
        if (ptrace(PT_CONTINUE, pid, 0, 0) != 0)
        {
          printf("tracer for %d ptrace continue without SIGSEGV\n", pid);
          tmp = -4;
        }
      }
      else
      {
        if (ptrace(PT_CONTINUE, pid, 0, 0) != 0)
        {
          printf("tracer for %d ptrace continue without SIGSEGV\n", pid);
          tmp = -4;
        }
      }
      //      break;
    }


  } while (tmp == 0);

  printf("tracer for %d Exit code: %d\n", pid, tmp);

  if (tmp != -6)
  {
    if (ptrace(PT_DETACH, pid, 0, 0) != 0)
      printf("tracer for %d ptrace detach", pid);
  }

  return 0;
}

void *tracer_thread(void *pid_ptr)
{
  isMain = 0;
  pid_t pid;
  pid = * ((pid_t *) pid_ptr);

  traceloop(pid);
  return 0;
}

void startThreadTracers(pid_t pid)
{
  FetchMap(pid);
  int j;
  for (j = 1 ; j < numberOfThreads; j++)
  {
    printf("tracer for %d start threads[%d] = %d\n", pid, j, threads[j]);
    assert(0 == pthread_create(&tids[j], NULL, tracer_thread, &threads[j]));
  }
}

void waitForAllTracers(pid_t pid)
{
  int j;
  for (j = 1 ; j < numberOfThreads; j++)
  {
    printf("tracer for %d join threads[%d] = %d\n", pid, j, threads[j]);
    assert(0 == pthread_join(tids[j], NULL));
    printf("tracer for %d after join threads[%d] = %d\n", pid, j, threads[j]);
  }
}

void openFilter()
{
  status = pipe(fildes);

  filter = fork();

  if (filter == 0) //child process
  {
    close(fildes[1]);
    unsigned long code;
    fd_set rfds;
    FD_ZERO(&rfds);
    FD_SET(fildes[0], &rfds);

    while (1)
    {
      int len;
      select(fildes[0] + 1, &rfds, NULL, NULL, NULL);
      if ((len = read(fildes[0], &code, sizeof(code))) < 0) break;
      else
      {
        if (len > 0)
        {
          if (code > 0) insertAddress(code);
	  else
            break;
	}
      }
      len = 0;
    }
    close(fildes[0]);
    outputAddresses(cmd);
    generategdbscripts();
    printf("tracer: filter stopping...\n");
    exit(0);
  }
  else
  {
    close(fildes[0]);
  }
}

void closeFilter()
{
  printf("tracer stopping...\n");
  unsigned long exitCode = 0;

  printf("tracer: notifying filter to exit\n");
  write(fildes[1], &exitCode, sizeof(exitCode));

  waitpid(filter, NULL, 0);
  printf("tracer: filter stopped...%d\n", filter);
  close(fildes[1]);
}

int main(int argc, char *argv[])
{
  isMain = 1;
  processParameters(argc, argv);
  openFilter();
  startThreadTracers(threads[0]);
  traceloop(threads[0]);
  waitForAllTracers(threads[0]);

  closeFilter();
  return 0;
}
