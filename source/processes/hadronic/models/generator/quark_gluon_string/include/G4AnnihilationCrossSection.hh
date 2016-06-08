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
