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
  public:
    G4BetaPlusDecayChannel (G4int Verbose,
                            const G4ParticleDefinition *theParentNucleus,
                            G4double theBR,
                            G4double theEndPointEnergy=0.0,
                            G4double theDaughterExcitation=0.0,
                            G4double theFFN=0.0) :
      G4NuclearDecayChannel (BetaPlus, Verbose, theParentNucleus, theBR, theFFN,
                             theEndPointEnergy,
                             theParentNucleus->GetBaryonNumber(),
                             int(theParentNucleus->GetPDGCharge()/eplus)-1,
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

