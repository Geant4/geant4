#ifndef G4AlphaDecayChannel_h
#define G4AlphaDecayChannel_h 1
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              G4AlphaDecayChannel.hh
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
class G4AlphaDecayChannel : public G4NuclearDecayChannel
{
  // class description 
  //
  //   Derived class from G4NuclearDecayChannel.  It is specific for
  //   Alpha decay proceess. 
  //
  // class  description - end
  public:
    G4AlphaDecayChannel (G4int Verbose,
                         const G4ParticleDefinition *theParentNucleus,
                         G4double theBR,
                         G4double theEndPointEnergy=0.0,
                         G4double theDaughterExcitation=0.0) :
      G4NuclearDecayChannel (Alpha, Verbose, theParentNucleus, theBR,
                             theEndPointEnergy,
                             (theParentNucleus->GetBaryonNumber())-4,
                             int(theParentNucleus->GetPDGCharge()/eplus)-2,
                             theDaughterExcitation, "alpha")
    {
#ifdef G4VERBOSE
      if (GetVerboseLevel()>1)
        G4cout <<"G4AlphaDecayChannel constructor" << G4endl;
#endif
    }
    ~G4AlphaDecayChannel () {;}
};
#endif

