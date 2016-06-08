#include "G4MesonSplitter.hh"

G4bool G4MesonSplitter::SplitMeson(G4int PDGcode, G4int* aEnd, G4int* bEnd)
{
  G4int absPDGcode = abs(PDGcode);
  if (absPDGcode >= 1000) return FALSE;

  G4int heavy =  absPDGcode/100;
  G4int light = (absPDGcode%100)/10;
  G4int anti  = 1 - 2*(max(heavy, light)%2);
  if (PDGcode < 0 ) anti = -anti;
  heavy *=  anti;
  light *= -anti;
  if ( anti < 0) 
     G4SwapObj(&heavy, &light);
  *aEnd = heavy;
  *bEnd = light;
  return TRUE; 	 
}
