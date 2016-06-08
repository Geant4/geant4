#include "G4ASCCrossSection.hh"

G4ASCCrossSection::
G4ASCCrossSection(G4int aCode1, G4int aCode2, G4double aX,  G4double aY, 
                  G4double aEta, G4double aEps)
{
  theCode1 = aCode1;
  theCode2 = aCode2;
  theX = aX;
  theY = aY;
  theEta = aEta;
  theEps = aEps;
}
