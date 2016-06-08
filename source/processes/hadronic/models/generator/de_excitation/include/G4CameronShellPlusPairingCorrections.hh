#ifndef G4CameronShellPlusPairingCorrections_h
#define G4CameronShellPlusPairingCorrections_h 1

#include "globals.hh"

class G4CameronShellPlusPairingCorrections
{
private:

	// Dummy constructor
	G4CameronShellPlusPairingCorrections(G4double dummy);

	static G4CameronShellPlusPairingCorrections theInstance;
	
public:
	
	~G4CameronShellPlusPairingCorrections() {};

	static G4double GetShellPlusPairingZ(const G4int Z) {
		if (Z <= TableSize && Z > 1) return SPZTable[Z-1]*MeV;
		else {
			G4cerr << "G4CameronShellPlusPairingCorrections: out of table for Z = " << Z << G4endl;
			return 0.0;
		}
	}
	
	static G4double GetShellPlusPairingN(const G4int N) {
		if (N <= TableSize && N > 0) return SPNTable[N-1]*MeV;
		else {
			G4cerr << "G4CameronShellPlusPairingCorrections: out of table for N = " << N << G4endl;
			return 0.0;
		}
	}
	
		
	enum  { TableSize = 200 };

private:

	static G4double SPZTable[TableSize];

	static G4double SPNTable[TableSize];
	
};
#endif
