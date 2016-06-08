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
// $Id: G4CameronGilbertShellCorrections.hh,v 1.5 2001/08/01 17:04:09 hpw Exp $
// GEANT4 tag $Name: geant4-04-00 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara

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
