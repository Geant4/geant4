#ifndef G4ASCCrossSection_h
#define G4ASCCrossSection_h

#include "globals.hh"
#include "G4VAnnihilationCrossSection.hh"

class G4ASCCrossSection : public G4VAnnihilationCrossSection
{
  public:
    G4ASCCrossSection(G4int, G4int, G4double,  G4double, G4double, G4double);
    G4bool InCharge(G4int aCode, G4int bCode);
    G4double GetXsec(G4double s);
  private:
    
    G4int theCode1;
    G4int theCode2;
    G4double theX;
    G4double theY;
    G4double theEta;
    G4double theEps;
};


inline G4bool G4ASCCrossSection::
InCharge(G4int aCode, G4int bCode)
{
  G4bool result;
  result = (aCode==theCode1&&bCode==theCode2)||(aCode==theCode2&&bCode==theCode1);
  return result;
}

inline G4double G4ASCCrossSection::
GetXsec(G4double s)
{
   G4double result = theX*pow(s, theEps) + theY*pow(s, -theEta);
   return result;
}

#endif
