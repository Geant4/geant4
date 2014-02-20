//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
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

G4ThreadLocal std::vector<sighandler_t> *G4HadSignalHandler::theCache = 0;
G4ThreadLocal bool G4HadSignalHandler::registered = false;

G4HadSignalHandler::G4HadSignalHandler(sighandler_t aNew)
{
    if (!theCache) theCache = new std::vector<sighandler_t>;
    if(!registered) 
    { 
      G4HadSignalHandler_local::G4HadSignalHandler_initial = 
         signal(SIGSEGV, G4HadSignalHandler_local::HandleIt);
      registered = true;
    }
    theCache->push_back(aNew);
}

G4HadSignalHandler::~G4HadSignalHandler()
{
  theCache->clear();
  signal (SIGSEGV, G4HadSignalHandler_local::G4HadSignalHandler_initial); 
  registered = false;
}

void G4HadSignalHandler_local::HandleIt(int i)
{
  static G4ThreadLocal int *iii_p = 0 ;
  if (!iii_p)
  {
    iii_p = new int ;
    *iii_p = G4HadSignalHandler::theCache->size()-1 ;
  }
  int &iii = *iii_p;
  for(int c=iii; c!=-1; c--)
  {
    iii--;
    //Andrea Dotti (13Jan2013): change for G4MT
    (G4HadSignalHandler::theCache->operator[](c))(i);
    //G4HadSignalHandler::theCache[c](i);
  }
    std::cerr << "callback to user-defined or default signal handler"<<endl;
  signal (SIGSEGV, G4HadSignalHandler_local::G4HadSignalHandler_initial);
  raise(i);
}

#endif
