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
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
#ifndef G4BetaMinusDecayChannel_h
#define G4BetaMinusDecayChannel_h 1
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              G4BetaMinusDecayChannel.hh
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
#include "globals.hh"
#include "G4NuclearDecayChannel.hh"
#include "G4RadioactiveDecayMode.hh"
  ////////////////////////////////////////////////////////////////////////////////
//
class G4BetaMinusDecayChannel : public G4NuclearDecayChannel 
{
  // class description 
  //
  //   Derived class from G4NuclearDecayChannel.  It is specific for
  //   Beta- decay proceess. 
  //
  // class  description - end
  public:
    G4BetaMinusDecayChannel (G4int Verbose,
                             const G4ParticleDefinition *theParentNucleus,
                             G4double theBR,
                             G4double theEndPointEnergy=0.0,
                             G4double theDaughterExcitation=0.0,
                             G4double theFFN=1.0,
			     G4bool   theBetaSimple = false,
			     RandGeneral* theRandEnergy = NULL):
      G4NuclearDecayChannel (BetaMinus, Verbose, theParentNucleus, theBR,
                             theFFN, 
			     theBetaSimple,
			     theRandEnergy,
			     theEndPointEnergy,
                             theParentNucleus->GetBaryonNumber(),
                             int(theParentNucleus->GetPDGCharge()/eplus)+1,
                             theDaughterExcitation, "e-", "anti_nu_e")
    {
#ifdef G4VERBOSE
      if (GetVerboseLevel()>1)
        G4cout <<"G4BetaMinusDecayChannel constructor" <<G4endl;
#endif
    }
    ~G4BetaMinusDecayChannel () {;}
};
#endif

