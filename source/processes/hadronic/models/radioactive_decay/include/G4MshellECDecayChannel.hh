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
#ifndef G4MshellECDecayChannel_h
#define G4MshellECDecayChannel_h 1

#include "globals.hh"
#include "G4NuclearDecayChannel.hh"
#include "G4RadioactiveDecayMode.hh"
////////////////////////////////////////////////////////////////////////////////
//
class G4MshellECDecayChannel : public G4NuclearDecayChannel
{
  // class description 
  //
  //   Derived class from G4NuclearDecayChannel.  It is specific for
  //   Electron Capture proceess with electron from the M shell
  //
  // class  description - end

public:
    G4MshellECDecayChannel (G4int Verbose,
			    const G4ParticleDefinition *theParentNucleus, 
			    G4double theBR,G4double theQtransit=0.0,
			    G4double theDaughterExcitation=0.0) :
      G4NuclearDecayChannel (MshellEC, Verbose, theParentNucleus, theBR,
			     theQtransit, theParentNucleus->GetBaryonNumber(),
			     int(theParentNucleus->GetPDGCharge()/eplus)-1, 
			     theDaughterExcitation,
			     "nu_e")
  {  
#ifdef G4VERBOSE
    if (GetVerboseLevel()>1)
      G4cerr <<"G4MshellECDecayChannel constructor" << G4endl;
#endif
  }
  ~G4MshellECDecayChannel () {;}
};
#endif

