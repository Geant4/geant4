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
#ifndef G4HadSignalHandler_hh
#define G4HadSignalHandler_hh

#ifndef G4HadSignalHandler_off

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

#endif
