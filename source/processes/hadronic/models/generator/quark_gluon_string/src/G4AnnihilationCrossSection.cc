#include "G4AnnihilationCrossSection.hh"
#include "G4ASCCrossSection.hh"
#include "G4GammaAnnCrossSection.hh"

G4AnnihilationCrossSection::
G4AnnihilationCrossSection()
{
// pi+ p
  theDataSets.push_back(new G4ASCCrossSection(211,  2212, 13.7, 27.8, 0.45, 0.079));
// pi+ n
  theDataSets.push_back(new G4ASCCrossSection(211,  2112, 13.7, 27.8, 0.45, 0.079));
// pi- p
  theDataSets.push_back(new G4ASCCrossSection(-211,  2212, 13.7, 35.9, 0.45, 0.079));
// pi- n
  theDataSets.push_back(new G4ASCCrossSection(-211,  2112, 13.7, 35.9, 0.45, 0.079));
// pi0 p
  theDataSets.push_back(new G4ASCCrossSection(111,  2212, 13.7, 35.9, 0.45, 0.079));
// pi0 n
  theDataSets.push_back(new G4ASCCrossSection(111,  2112, 13.7, 35.9, 0.45, 0.079));
// K-  p 
  theDataSets.push_back(new G4ASCCrossSection(-321,  2212, 12.2, 26.4, 0.50, 0.079));
// K-  n 
  theDataSets.push_back(new G4ASCCrossSection(-321,  2112, 12.2, 26.4, 0.50, 0.079));
// K0  p 
  theDataSets.push_back(new G4ASCCrossSection(-311,  2212, 12.2, 26.4, 0.50, 0.079));
// K0  n 
  theDataSets.push_back(new G4ASCCrossSection(-311,  2112, 12.2, 26.4, 0.50, 0.079));
// p-  p
  theDataSets.push_back(new G4ASCCrossSection(-2212,  2212, 22.0, 98.2, 0.46, 0.079));
// p-  n
  theDataSets.push_back(new G4ASCCrossSection(-2212,  2112, 22.3, 92.7, 0.46, 0.079));
// n-  n
  theDataSets.push_back(new G4ASCCrossSection(-2112,  2112, 22.0, 98.2, 0.46, 0.079));
// n-  p
  theDataSets.push_back(new G4ASCCrossSection(-2112,  2212, 22.3, 92.7, 0.46, 0.079));
// gamma+n and gamma+p
  theDataSets.push_back(new G4GammaAnnCrossSection);
}
