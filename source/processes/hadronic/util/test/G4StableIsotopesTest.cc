// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4StableIsotopesTest.cc,v 1.2 1999-12-15 14:53:41 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "G4ios.hh"
#include "g4std/fstream"
#include "g4std/iomanip"
#include "G4StableIsotopes.hh"

main()
{
   G4StableIsotopes theIso;
   for (int Z=1; Z<92; Z++)
   {
     G4int protonCount = theIso.GetProtonCount(Z);
     G4String theName = theIso.GetName(Z);
     for (G4int i=0; i<theIso.GetNumberOfIsotopes(Z); i++)
     {
       G4int nucleons = theIso.GetIsotopeNucleonCount(theIso.GetFirstIsotope(Z)+i);
       G4double fracInPercent=theIso.GetAbundance(theIso.GetFirstIsotope(Z)+i);
       G4cout << nucleons << " " << fracInPercent << " " << protonCount << " " << theName << G4endl;
     }
     G4cout << G4endl;
   }
   G4cout << G4endl;
}
