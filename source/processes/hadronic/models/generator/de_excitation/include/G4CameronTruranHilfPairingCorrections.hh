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
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4CameronTruranHilfPairingCorrections.hh,v 1.5 2001/08/01 17:04:10 hpw Exp $
// GEANT4 tag $Name: geant4-04-00 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara

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
