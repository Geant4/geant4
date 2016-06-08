#ifndef G4CameronGilbertPairingCorrections_h
#define G4CameronGilbertPairingCorrections_h 1

#include "globals.hh"

class G4CameronGilbertPairingCorrections
{
private:
	// Dummy constructor
	G4CameronGilbertPairingCorrections(G4double dummy);
	
	static G4CameronGilbertPairingCorrections theInstance;


public:
	
	~G4CameronGilbertPairingCorrections() {};

	G4double GetPairingZ(const G4int Z) {
		if (Z <= ZTableSize && Z > 1) return PairingZTable[Z-1]*MeV;
		else {
			G4cerr << "G4CameronGilbertPairingCorrections: out of table for Z = " << Z << G4endl;
			return 0.0;
		}
	}
	
	G4double GetPairingN(const G4int N) {
		if (N <= NTableSize && N > 0) return PairingNTable[N-1]*MeV;
		else {
			G4cerr << "G4CameronGilbertPairingCorrections: out of table for N = " << N << G4endl;
			return 0.0;
		}
	}
	
	
	enum  { ZTableSize = 98, NTableSize = 150 };
private:



	static G4double PairingZTable[ZTableSize];

	static G4double PairingNTable[NTableSize];
	
};
#endif
