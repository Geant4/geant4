
#include "globals.hh"


#ifndef G4NucleiPropertiesTheoreticalTable_h
#define G4NucleiPropertiesTheoreticalTable_h 1


//	Encapsulates Data from W.D. Myers, W.J. Swiatecki, P. Moller and J.R. Nix, 
//	1. Jan. 1995.
//	Atomic Mass Excess
// 


class G4NucleiPropertiesTheoreticalTable 
{
private:
  
  // Default constructor 
  G4NucleiPropertiesTheoreticalTable(G4double dummy);

  static G4NucleiPropertiesTheoreticalTable theInstance;

public:

  // Destructor
  ~G4NucleiPropertiesTheoreticalTable() { };

  enum  {nEntries = 8979, shortTableSize = 137}; 

  // Other Operations 

  // Operation: GetMassExcess
  static G4double GetMassExcess(G4int Z, G4int A); 


	// Operation: GetNuclearMass
	static G4double GetNuclearMass(G4int Z, G4int A);


	// Operation: GetAtomicMass 
	static G4double GetAtomicMass(G4int Z, G4int A);

	// Operation: GetBindingEnergy
	static G4double GetBindingEnergy(G4int Z, G4int A);


	// Is the nucleus (Z,A) in table?
	static G4bool IsInTable(G4int Z, G4int A);



private:

	// Operation: GetIndex
	static G4int GetIndex(G4int Z, G4int A);
  
	static G4double ElectronicBindingEnergy(G4int Z);
 



	// Mass Excess
	static G4double AtomicMassExcess[nEntries];
  
  

    
	// Table of Z (number of protons) and A (number of nucleons)
	//        indexArray[0][ ] --> Z
	//        indexArray[1][ ] --> A
	static G4int indexArray[2][nEntries];

	// Reduced Table of Z for shorter index search.
	//         The index in this table coincide with Z-1
	//         For each Z value shortTable[Z-1] has the index of the 1st occurrence in
	//         the indexArray[][]
	static G4int shortTable[shortTableSize];


};
  
inline G4double G4NucleiPropertiesTheoreticalTable::GetMassExcess(G4int Z, G4int A) 
{
    G4int i=GetIndex(Z, A);
    if (i >= 0) {
		return AtomicMassExcess[i]*MeV;
    } else {
        return 0.0;
    }
}

inline G4double G4NucleiPropertiesTheoreticalTable::GetBindingEnergy(G4int Z, G4int A)
{
    G4int i=GetIndex(Z, A);
    if (i >= 0){
	    const G4double Mh = 7.289034*MeV;  // hydrogen atom mass excess
		 const G4double Mn = 8.071431*MeV;  // neutron mass excess
		 return G4double(Z)*Mh + G4double(A-Z)*Mn - AtomicMassExcess[i]*MeV;
    } else { 
	    return 0.0;
    }
}



inline G4double  G4NucleiPropertiesTheoreticalTable::GetAtomicMass(G4int Z, G4int A)
{
    G4int i=GetIndex(Z, A);
    if (i >= 0) {
      return AtomicMassExcess[i]*MeV + A*amu_c2;
    } else {
      return 0.0;
    }
}
  


inline G4double  G4NucleiPropertiesTheoreticalTable::GetNuclearMass(G4int Z, G4int A)
{
    G4int i=GetIndex(Z, A);
    if (i >= 0) {
      return GetAtomicMass(Z,A) - G4double(Z)*electron_mass_c2 + ElectronicBindingEnergy(Z);
    } else {
      return 0.0;
    }
}

inline G4double G4NucleiPropertiesTheoreticalTable::ElectronicBindingEnergy(G4int Z) {
	const G4double ael = 1.433e-5*MeV; // electronic-binding constant
	return ael*pow(G4double(Z),2.39);
}

inline G4bool G4NucleiPropertiesTheoreticalTable::IsInTable(G4int Z, G4int A)
{
    return (Z <= A && A >= 16 && A <= 339 && Z <= 136 && Z >= 8 && GetIndex(Z, A) >= 0);
}


#endif






