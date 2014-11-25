#ifndef G4DNAChemistryList_hh
#define G4DNAChemistryList_hh 1

#include "G4VPhysicsConstructor.hh"
#include "G4VUserChemistryList.hh"
#include "globals.hh"

class G4DNAMolecularReactionTable;

class G4EmDNAChemistry : public G4VUserChemistryList,
                         public G4VPhysicsConstructor
{

public:
  G4EmDNAChemistry();
  virtual ~G4EmDNAChemistry();

  virtual void ConstructParticle()
  {
    ConstructMolecule();
  }
  virtual void ConstructMolecule();
  virtual void ConstructProcess();

  virtual void ConstructDissociationChannels();
  virtual void ConstructReactionTable(G4DNAMolecularReactionTable* reactionTable);
  virtual void ConstructTimeStepModel(G4DNAMolecularReactionTable* reactionTable);

};

#endif
