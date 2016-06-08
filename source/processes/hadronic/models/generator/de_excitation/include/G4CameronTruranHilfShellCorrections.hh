#ifndef G4CameronTruranHilfShellCorrections_h
#define G4CameronTruranHilfShellCorrections_h 1

#include "globals.hh"

class G4CameronTruranHilfShellCorrections
{
private:
	
	// Dummy constructor
	G4CameronTruranHilfShellCorrections(G4double dummy);
	
	static G4CameronTruranHilfShellCorrections theInstance;


public:
	
	~G4CameronTruranHilfShellCorrections() {};

	static G4double GetShellZ(const G4int Z) {
		if (Z <= ZTableSize && Z > 1) return ShellZTable[Z-1]*MeV;
		else {
			G4cerr << "G4CameronTruranHilfShellCorrections: out of table for Z = " << Z << G4endl;
			return 0.0;
		}
	}
	
	static G4double GetShellN(const G4int N) {
		if (N <= NTableSize && N > 0) return ShellNTable[N-1]*MeV;
		else {
			G4cerr << "G4CameronTruranHilfShellCorrections: out of table for N = " << N << G4endl;
			return 0.0;
		}
	}
	
	
	enum  { ZTableSize = 102, NTableSize = 155 };
private:



	static G4double ShellZTable[ZTableSize];

	static G4double ShellNTable[NTableSize];
	
};
#endif
