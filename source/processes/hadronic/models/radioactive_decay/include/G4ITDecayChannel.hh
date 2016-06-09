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
#ifndef G4ITDecayChannel_h
#define G4ITDecayChannel_h 1

#include "globals.hh"
#include "G4NuclearDecayChannel.hh"
////////////////////////////////////////////////////////////////////////////////
//
class G4ITDecayChannel : public G4NuclearDecayChannel 
{
  
  // class description 
  //
  //   Derived class from G4NuclearDecayChannel.  It is specific for
  //   Isomeric Transitions 
  //
  // class  description - end
  public:
    G4ITDecayChannel (G4int Verbose,
                      const G4Ions *theParentNucleus,
                      G4double theBR) :
      G4NuclearDecayChannel (IT, Verbose, theParentNucleus, theBR, 0.0,
			     theParentNucleus->GetBaryonNumber(),
			     int(theParentNucleus->GetPDGCharge()/eplus),
			     theParentNucleus->GetExcitationEnergy())
    {
#ifdef G4VERBOSE
      if (GetVerboseLevel()>1)
         G4cout <<"G4ITDecayChannel constructor" << G4endl;
#endif
    }
    ~G4ITDecayChannel () {;}
};
#endif

