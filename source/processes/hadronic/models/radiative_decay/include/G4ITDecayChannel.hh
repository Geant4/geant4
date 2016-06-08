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
                      const G4ParticleDefinition *theParentNucleus,
                      G4double theBR) :
      G4NuclearDecayChannel (IT, Verbose, theParentNucleus, theBR, 0.0,
			     theParentNucleus->GetBaryonNumber(),
			     int(theParentNucleus->GetPDGCharge()/eplus),
			     ((const G4Ions*) theParentNucleus)->GetExcitationEnergy())
    {
#ifdef G4VERBOSE
      if (GetVerboseLevel()>1)
         G4cout <<"G4ITDecayChannel constructor" << G4endl;
#endif
    }
    ~G4ITDecayChannel () {;}
};
#endif

