//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
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
