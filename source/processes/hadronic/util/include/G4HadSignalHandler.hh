#ifndef G4HadSignalHandler_hh
#define G4HadSignalHandler_hh

#include <iostream>
#include <signal.h>
#include <vector>

// A simple, reasonably portable, but at present 
// semantic-wise totally unsafe signalhandler prototype meant for SEGFAULT.
// Being rushed into production or various reasons.

extern "C"
{
  typedef void (*sighandler_t)(int);
}

class G4HadSignalHandler
{
  public: 
  
  G4HadSignalHandler(sighandler_t aNew);
  
  ~G4HadSignalHandler();
      
  static std::vector<sighandler_t> theCache;  
  static bool registered;
};

 

#endif
