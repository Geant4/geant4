#include "outputList.hh"
#include <strstream>
#include <iomanip>
#include <fstream>

outputList::~outputList() 
{ 
  for (int i=0; i<N; i++)
    if ( erase[i] ) 
      delete files[i];
  delete [] files; 
}

outputList::outputList(int n) 
  : N(n),files(new std::ostream*[n]),erase(new bool[n]) 
{
  for (int i=0; i<N; i++) {
    files[i] = new std::ofstream("/dev/null");
    erase[i] = true;
  }
}

outputList::outputList(int n,char* file) 
  : N(n),files(new std::ostream*[n]),erase(new bool[n])
{
  int l = strlen(file)+4;
  char* s = new char[l];
  for (int i=0; i<N; i++) {
    if ( strcmp(file,"/dev/null") ) { 
      std::ostrstream string(s,l);
      string << file << "." << std::setw(2) << std::setfill('0') << i;
      files[i] = new std::ofstream(s);
    }
    else 
      files[i] = new std::ofstream("/dev/null");
    erase[i] = true;
  }
  delete s;
}

void outputList::add(int i,std::ostream* f)
{ 
  if (f) {
    if ( files[i] ) 
      delete files[i];
    files[i] = f;
    erase[i] = false;
  }
}

