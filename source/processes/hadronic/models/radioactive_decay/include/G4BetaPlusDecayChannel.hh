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
#ifndef G4BetaPlusDecayChannel_h
#define G4BetaPlusDecayChannel_h 1
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              G4BetaPlusDecayChannel.hh
//
// Version:             0.b.3
// Date:                29/02/00
// Author:              F Lei & P R Truscott
// Organisation:        DERA UK
// Customer:            ESA/ESTEC, NOORDWIJK
// Contract:            12115/96/JG/NL Work Order No. 3
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// CHANGE HISTORY
// --------------
//
// 29 February 2000, P R Truscott, DERA UK
// 0.b.3 release.
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
////////////////////////////////////////////////////////////////////////////////
//

#include <CLHEP/Units/SystemOfUnits.h>

#include "globals.hh"
#include "G4NuclearDecayChannel.hh"
#include "G4RadioactiveDecayMode.hh"

////////////////////////////////////////////////////////////////////////////////
//
class G4BetaPlusDecayChannel : public G4NuclearDecayChannel
{
  // class description 
  //
  //   Derived class from G4NuclearDecayChannel.  It is specific for
  //   Beta+ decay proceess. 
  //
  // class  description - end
  //
  //
  public:
    G4BetaPlusDecayChannel (G4int Verbose,
                            const G4ParticleDefinition *theParentNucleus,
                            G4double theBR,
                            G4double theEndPointEnergy=0.0,
                            G4double theDaughterExcitation=0.0,
                            G4double theFFN=0.0,
			    G4bool   theBetaSimple = false,
			    CLHEP::RandGeneral* theRandEnergy = 0) :
      G4NuclearDecayChannel (BetaPlus, 
			     Verbose, 
			     theParentNucleus, 
			     theBR,
			     theFFN, 
			     theBetaSimple,
			     theRandEnergy,
                             theEndPointEnergy,
                             theParentNucleus->GetBaryonNumber(),
                             int(theParentNucleus->GetPDGCharge()/CLHEP::eplus)-1,
                             theDaughterExcitation, "e+", "nu_e")
    {
#ifdef G4VERBOSE
      if (GetVerboseLevel()>1)
        G4cout <<"G4BetaPlusDecayChannel constructor" <<G4endl;
#endif
    }
  ~G4BetaPlusDecayChannel () {;}
};
#endif

