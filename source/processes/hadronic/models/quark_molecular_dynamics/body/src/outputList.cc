#include "outputList.hh"
#include <strstream.h>
#include <iomanip.h>
#include <fstream.h>

outputList::~outputList() 
{ 
  for (int i=0; i<N; i++)
    if ( erase[i] ) 
      delete files[i];
  delete [] files; 
}

outputList::outputList(int n) 
  : N(n),files(new ostream*[n]),erase(new bool[n]) 
{
  for (int i=0; i<N; i++) {
    files[i] = new ofstream("/dev/null");
    erase[i] = true;
  }
}

outputList::outputList(int n,char* file) 
  : N(n),files(new ostream*[n]),erase(new bool[n])
{
  int l = strlen(file)+4;
  char* s = new char[l];
  for (int i=0; i<N; i++) {
    if ( strcmp(file,"/dev/null") ) { 
      ostrstream string(s,l);
      string << file << "." << setw(2) << setfill('0') << i;
      files[i] = new ofstream(s);
    }
    else 
      files[i] = new ofstream("/dev/null");
    erase[i] = true;
  }
  delete s;
}

void outputList::add(int i,ostream* f)
{ 
  if (f) {
    if ( files[i] ) 
      delete files[i];
    files[i] = f;
    erase[i] = false;
  }
}

