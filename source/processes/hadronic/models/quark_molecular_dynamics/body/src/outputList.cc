#include "outputList.hh"
#include "g4std/strstream"
#include "g4std/iomanip"
#include "g4std/fstream"

outputList::~outputList() 
{ 
  for (int i=0; i<N; i++)
    if ( erase[i] ) 
      delete files[i];
  delete [] files; 
}

outputList::outputList(int n) 
  : N(n),files(new G4std::ostream*[n]),erase(new bool[n]) 
{
  for (int i=0; i<N; i++) {
    files[i] = new G4std::ofstream("/dev/null");
    erase[i] = true;
  }
}

outputList::outputList(int n,char* file) 
  : N(n),files(new G4std::ostream*[n]),erase(new bool[n])
{
  int l = strlen(file)+4;
  char* s = new char[l];
  for (int i=0; i<N; i++) {
    if ( strcmp(file,"/dev/null") ) { 
      G4std::ostrstream string(s,l);
      string << file << "." << G4std::setw(2) << setfill('0') << i;
      files[i] = new G4std::ofstream(s);
    }
    else 
      files[i] = new G4std::ofstream("/dev/null");
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

