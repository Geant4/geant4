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
#ifndef G4AnnihilationCrossSection_h
#define G4AnnihilationCrossSection_h

#include "globals.hh"
#include "g4std/vector"  
#include "G4VAnnihilationCrossSection.hh"

class G4AnnihilationCrossSection
{
  public:
    G4AnnihilationCrossSection();
    
    G4double GetCrossSection(int aCode, int bCode, G4double s);
  private:
  
  G4std::vector<G4VAnnihilationCrossSection*> theDataSets;
};

inline G4double G4AnnihilationCrossSection::GetCrossSection(G4int aCode, G4int bCode, G4double s)
{
  G4double result = 0.;
  typedef G4std::vector<G4VAnnihilationCrossSection*>::iterator iter;
  iter i;
  for(i=theDataSets.begin(); i!=theDataSets.end(); i++)
  {
    if((*i)->InCharge(aCode, bCode))
    {
      result = (*i)->GetXsec(s);
      break;
    }
  }
  return result;
}
#endif
