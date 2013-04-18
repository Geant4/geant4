#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <unistd.h>
#include "memoryprotection/filter.h"

int position = -1;

//#define ADDRESSTABLEDEBUG 1

struct AddressTable
{
  int counter;
  int space;
  unsigned long *addresses;
} Detected;

int bSearchAddress(unsigned long address, int head, int tail)
{
  if (head >= tail) 
  {
    position = tail;
      
    return -1;
  }  
  int middle = (head + tail) / 2;
  if (Detected.addresses[middle] == address)
  {
    position = middle;
    return 0;
  }
  else if (Detected.addresses[middle] < address)
    return bSearchAddress(address, middle + 1, tail);
  else
    return bSearchAddress(address, 0, middle);
}

int findAddress(unsigned long address)
{
  if (position == -1) return -1;
  return bSearchAddress(address, 0, Detected.counter);
}

int insertAddress(unsigned long address)
{
  if (position == -1)
  {
    Detected.space = 4096;
    assert(0 != (Detected.addresses = (unsigned long *) malloc(sizeof(unsigned long) * Detected.space)));
    Detected.counter = 1;
    position = 0;
    #ifdef ADDRESSTABLEDEBUG
    printf("Inserted position: %d, address: 0x%lx.\n", position, address);
    #endif
    Detected.addresses[position] = address;
    return 0;
  }
  else
  {
    if (-1 == findAddress(address))
    {
      #ifdef ADDRESSTABLEDEBUG
      printf("Inserted position: %d, address: 0x%lx.\n", position, address);
      #endif
      Detected.counter++;
      if (Detected.counter > Detected.space)
      {
        Detected.space = Detected.space + 4096;

        assert(0 != (Detected.addresses = (unsigned long *) realloc((char *)Detected.addresses, sizeof(unsigned long) * Detected.space)));
      }
      int i;
      for (i = Detected.counter - 1; i > position ; i--)
        Detected.addresses[i] = Detected.addresses[i - 1];
      Detected.addresses[i] = address;

      return 0;
    }
    else
      return -1;
  }
}

void outputAddresses(char *cmd)
{
  int i;
  char convertcmd[1024];
  printf("Lines that produce segmentation faults:\n");
  for (i = 0 ; i < Detected.counter ; i++)
  {
    sprintf(convertcmd, "addr2line -e %s 0x%lx", cmd, Detected.addresses[i]);
    system(convertcmd);
    printf("%s address[%d]= 0x%lx\n", cmd, i, Detected.addresses[i]); 
  }
}

void generategdbscripts()
{
  int i;

  pid_t myselfpid = getpid();
  char gdbcommands[1024];
  gdbcommands[0]  = 0;
  sprintf(gdbcommands, "./gdbcommands.%d", myselfpid);

  FILE *fp = fopen(gdbcommands, "w+");
  for (i = 0 ; i < Detected.counter ; i++)
  {
    fprintf(fp, "break *0x%lx\n", Detected.addresses[i]); 
    fprintf(fp, "commands\n");
    fprintf(fp, "bt 10\n");
    fprintf(fp, "cont\n");
    fprintf(fp, "end\n");
  }
  fprintf(fp, "cont\n");
  fclose(fp);
}

#ifdef ADDRESSTABLEDEBUG
int main(int argc, char *argv[])
{
  insertAddress(1423423);
  insertAddress(1238743);
  insertAddress(2342134);
  outputAddresses();  
}
#endif
