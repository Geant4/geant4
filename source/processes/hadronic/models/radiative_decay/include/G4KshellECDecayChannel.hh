#ifndef G4KshellECDecayChannel_h
#define G4KshellECDecayChannel_h 1
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              G4KshellECDecayChannel.hh
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
class G4KshellECDecayChannel : public G4NuclearDecayChannel
{
  // class description 
  //
  //   Derived class from G4NuclearDecayChannel.  It is specific for
  //   Electron Capture proceess with electron from the K shell
  //
  // class  description - end

  public:
    G4KshellECDecayChannel (G4int Verbose,
                            const G4ParticleDefinition *theParentNucleus,
                            G4double theBR,
                            G4double theQtransit=0.0,
                            G4double theDaughterExcitation=0.0) :
      G4NuclearDecayChannel (KshellEC, Verbose, theParentNucleus, theBR,
                             theQtransit, theParentNucleus->GetBaryonNumber(),
                             int(theParentNucleus->GetPDGCharge()/eplus)-1,
                             theDaughterExcitation, "nu_e")
    {
#ifdef G4VERBOSE
      if (GetVerboseLevel()>1)
        G4cout <<"G4KshellECDecayChannel constructor" << G4endl;
#endif
    }
  ~G4KshellECDecayChannel () {;}
};
#endif

