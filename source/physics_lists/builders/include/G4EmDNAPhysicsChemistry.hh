#ifndef G4DNAChemistryList_hh
#define G4DNAChemistryList_hh 1

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"

class G4DNAMolecularReactionTable;

class G4EmDNAPhysicsChemistry :  public G4VPhysicsConstructor
{

public :
    G4EmDNAPhysicsChemistry(G4int ver = 1);
    virtual ~G4EmDNAPhysicsChemistry();

    virtual void ConstructParticle();
    virtual void ConstructProcess();

protected :
    void ConstructMolecules() ;
    void ConstructDecayChannels() ;
    void ConstructReactionTable() ;

private:
  G4int  verbose;
};

#endif
