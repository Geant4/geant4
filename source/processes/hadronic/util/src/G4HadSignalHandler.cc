#ifndef G4HadSignalHandler_off

#include "G4HadSignalHandler.hh"


namespace G4HadSignalHandler_local
{
  extern "C" 
  {
    void HandleIt(int i);
    static void (*G4HadSignalHandler_initial)(int);
  }
}

using namespace std;

std::vector<sighandler_t> G4HadSignalHandler::theCache;
bool G4HadSignalHandler::registered = false;

G4HadSignalHandler::G4HadSignalHandler(sighandler_t aNew)
{
    if(!registered) 
    { 
      G4HadSignalHandler_local::G4HadSignalHandler_initial = 
         signal(SIGSEGV, G4HadSignalHandler_local::HandleIt);
      registered = true;
    }
    theCache.push_back(aNew);
}

G4HadSignalHandler::~G4HadSignalHandler()
{
  theCache.clear();
  signal (SIGSEGV, G4HadSignalHandler_local::G4HadSignalHandler_initial); 
  registered = false;
}

void G4HadSignalHandler_local::HandleIt(int i)
{
  static int iii=G4HadSignalHandler::theCache.size()-1;
  for(int c=iii; c!=-1; c--)
  {
    iii--;
    G4HadSignalHandler::theCache[c](i);
  }
    std::cerr << "callback to user-defined or default signal handler"<<endl;
  signal (SIGSEGV, G4HadSignalHandler_local::G4HadSignalHandler_initial);
  raise(i);
}

#endif
