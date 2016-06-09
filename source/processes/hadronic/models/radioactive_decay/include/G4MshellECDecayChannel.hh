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

