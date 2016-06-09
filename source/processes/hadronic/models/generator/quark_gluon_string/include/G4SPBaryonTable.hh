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
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
#ifndef G4SPBaryonTable_h
#define G4SPBaryonTable_h

#include <vector>
#include "G4SPBaryon.hh"

class G4SPBaryonTable
{
  public:
  struct DeleteSPBaryon{void operator()(G4SPBaryon* aS){delete aS;} };
  
  ~G4SPBaryonTable() {std::for_each(theBaryons.begin(), theBaryons.end(), G4SPBaryonTable::DeleteSPBaryon());}
  void insert(G4SPBaryon * aBaryon) { theBaryons.push_back(aBaryon);}
  G4double length() {return theBaryons.size();}
  
  const G4SPBaryon * GetBaryon(G4ParticleDefinition * aDefinition);

  private:
  std::vector<G4SPBaryon *> theBaryons;

};

inline const G4SPBaryon * G4SPBaryonTable::
GetBaryon(G4ParticleDefinition * aDefinition)
{
  G4SPBaryon * result = 0;
  for(unsigned int i=0; i<theBaryons.size(); i++)
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
