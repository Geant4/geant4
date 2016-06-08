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
// $Id: G4CameronTruranHilfShellCorrections.hh,v 1.3.2.1 2001/06/28 19:12:59 gunter Exp $
// GEANT4 tag $Name:  $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara

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
