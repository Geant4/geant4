//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4CameronGilbertPairingCorrections.hh,v 1.3.2.1 2001/06/28 19:12:58 gunter Exp $
// GEANT4 tag $Name:  $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara


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
