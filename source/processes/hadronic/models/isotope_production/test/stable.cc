#include "G4StableIsotopes.hh"
#include "globals.hh"
#include "g4std/iostream"

int main()
{
// To loop over all isotopes for element with protonount Z
 G4int Z;
 G4StableIsotopes theIso;
 for(Z=1; Z<93; Z++)
 {
   for (G4int i=0; i<theIso.GetNumberOfIsotopes(Z); i++)
   {
      cout <<theIso.GetName(Z)<<" ";
      cout <<Z<<" ";
      cout <<theIso.GetIsotopeNucleonCount(theIso.GetFirstIsotope(Z)+i)<<" ";
      cout << G4endl;
    }
    cout << G4endl;
  }
}
