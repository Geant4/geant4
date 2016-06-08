#ifndef G4CameronGilbertShellCorrections_h
#define G4CameronGilbertShellCorrections_h 1

#include "globals.hh"

class G4CameronGilbertShellCorrections
{
private:
	// Dummy constructor
	G4CameronGilbertShellCorrections(G4double dummy);
	
	static G4CameronGilbertShellCorrections theInstance;


public:
	
	~G4CameronGilbertShellCorrections() {};

	G4double GetShellZ(const G4int Z) const {
		if (Z <= ZTableSize && Z > 1) return ShellZTable[Z-1]*MeV;
		else {
			G4cerr << "G4CameronGilbertShellCorrections: out of table for Z = " << Z << G4endl;
			return 0.0;
		}
	}
	
	G4double GetShellN(const G4int N) const {
		if (N <= NTableSize && N > 0) return ShellNTable[N-1]*MeV;
		else {
			G4cerr << "G4CameronGilbertShellCorrections: out of table for N = " << N << G4endl;
			return 0.0;
		}
	}
	
	
	enum  { ZTableSize = 98, NTableSize = 150 };
private:



	static const G4double ShellZTable[ZTableSize];

	static const G4double ShellNTable[NTableSize];
	
};
#endif
