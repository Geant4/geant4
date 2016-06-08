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
// $Id: G4PairingCorrection.hh,v 1.1.2.1 2001/06/28 19:13:03 gunter Exp $
// GEANT4 tag $Name:  $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara

#ifndef G4PairingCorrection_h
#define G4PairingCorrection_h 1

#include "globals.hh"

class G4PairingCorrection
{
private:

  // Dummy constructor
  G4PairingCorrection(G4double dummy);

  static G4PairingCorrection theInstance;
   
public:
	
  ~G4PairingCorrection() {};

  static G4double GetPairingCorrection(const G4int anA, const G4int aZ) 
  {
    const G4double PairingConstant = 12.0*MeV;
    const G4int N = anA - aZ;
    G4double Pair = (1.0 - G4double(aZ) + 2.0*(aZ/2)) + (1.0 - G4double(N) + 2.0*(N/2));
    G4double PCorrection = Pair*PairingConstant/sqrt(G4double(anA));
    return PCorrection;
  }


  static G4double GetFissionPairingCorrection(const G4int A, const G4int Z)
  {
    const G4double PairingConstant = 14.0*MeV;
    const G4int N = A - Z;
    G4double Pair = (1.0 - G4double(Z) + 2.0*(Z/2)) + (1.0 - G4double(N) + 2.0*(N/2));
    G4double PCorrection = Pair*PairingConstant/sqrt(G4double(A));
    return PCorrection;
  }

	
};
#endif
