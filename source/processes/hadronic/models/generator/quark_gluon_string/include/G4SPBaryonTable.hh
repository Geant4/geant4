#ifndef G4SPBaryonTable_h
#define G4SPBaryonTable_h

#include "g4rw/tpordvec.h"
#include "G4SPBaryon.hh"

class G4SPBaryonTable
{
  public:
  ~G4SPBaryonTable() {theBaryons.clearAndDestroy();}
  void insert(G4SPBaryon * aBaryon) { theBaryons.insert(aBaryon);}
  G4double length() {return theBaryons.length();}
  
  const G4SPBaryon * GetBaryon(G4ParticleDefinition * aDefinition);

  private:
  G4RWTPtrOrderedVector<G4SPBaryon> theBaryons;

};

inline const G4SPBaryon * G4SPBaryonTable::
GetBaryon(G4ParticleDefinition * aDefinition)
{
  G4SPBaryon * result = 0;
  for(G4int i=0; i<theBaryons.length(); i++)
  {
    if(theBaryons[i]->GetDefinition()==aDefinition)
    {
      result = theBaryons[i];
      break;
    }
  }
  return result;
}

#endif
