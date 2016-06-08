
#include "G4NucleiPropertiesTheoreticalTable.hh"



G4int G4NucleiPropertiesTheoreticalTable::shortTable[G4NucleiPropertiesTheoreticalTable::shortTableSize];

G4NucleiPropertiesTheoreticalTable G4NucleiPropertiesTheoreticalTable::theInstance(2.5);
// Default constructor
G4NucleiPropertiesTheoreticalTable::G4NucleiPropertiesTheoreticalTable(G4double dummy)
{
	G4double even_more_dummy = dummy;
	G4int j = 0;
	
	for (G4int i = 0; i < G4NucleiPropertiesTheoreticalTable::nEntries; i++) {
		if (indexArray[0][i] > j+7) {
			shortTable[j++] = i;
		}
	}
	shortTable[G4NucleiPropertiesTheoreticalTable::shortTableSize-1] = 
										G4NucleiPropertiesTheoreticalTable::nEntries;
}

// Determine the table index for a Nuclide with Z protons and A nucleons
G4int G4NucleiPropertiesTheoreticalTable::GetIndex(G4int Z, G4int A) 
{

	if(A>339) G4Exception(
	"G4NucleiPropertiesTheoreticalTable::GetIndex : Nucleon number larger than 339!\n"); 
	if(A<16) G4Exception(
	"G4NucleiPropertiesTheoreticalTable::GetIndex : Nucleon number smaller than 16!\n"); 
	if(Z>136) G4Exception(
	"G4NucleiPropertiesTheoreticalTable::GetIndex : Proton number larger than 136!\n"); 
	if(Z<8) G4Exception(
	"G4NucleiPropertiesTheoreticalTable::GetIndex : Proton number smaller than 8!\n");	
	if(Z>A) G4Exception(
	"G4NucleiPropertiesTheoreticalTable::GetIndex : Nucleon number smaller than Z!\n"); 
  
  
  G4int i = shortTable[Z-8];
  while ( i < shortTable[Z-8+1] ) {
    if (indexArray[1][i] != A ) i++;
    else return i;
  }
  
  return -1;
}



