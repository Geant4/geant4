#ifndef G4CameronTruranHilfPairingCorrections_h
#define G4CameronTruranHilfPairingCorrections_h 1

#include "globals.hh"

class G4CameronTruranHilfPairingCorrections
{
private:
	// Dummy constructor
	G4CameronTruranHilfPairingCorrections(G4double dummy);
	
	static G4CameronTruranHilfPairingCorrections theInstance;


public:
	
	~G4CameronTruranHilfPairingCorrections() {};

	G4double GetPairingZ(const G4int Z) const {
		if (Z <= ZTableSize && Z > 1) return PairingZTable[Z-1]*MeV;
		else {
			G4cerr << "G4CameronTruranHilfPairingCorrections: out of table for Z = " << Z << G4endl;
			return 0.0;
		}
	}
	
	G4double GetPairingN(const G4int N) const {
		if (N <= NTableSize && N > 0) return PairingNTable[N-1]*MeV;
		else {
			G4cerr << "G4CameronTruranHilfPairingCorrections: out of table for N = " << N << G4endl;
			return 0.0;
		}
	}
	
	
	enum  { ZTableSize = 102, NTableSize = 155 };
private:



	static const G4double PairingZTable[ZTableSize];

	static const G4double PairingNTable[NTableSize];
	
};
#endif
