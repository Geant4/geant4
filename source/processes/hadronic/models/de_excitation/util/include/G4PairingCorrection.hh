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
// $Id: G4PairingCorrection.hh,v 1.1 2003/08/26 18:50:10 lara Exp $
// GEANT4 tag $Name: geant4-06-00-patch-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara

#ifndef G4PairingCorrection_h
#define G4PairingCorrection_h 1

#include "globals.hh"
#include "G4CookPairingCorrections.hh"
#include "G4CameronGilbertPairingCorrections.hh"
//#include "G4CameronTruranHilfPairingCorrections.hh"


class G4PairingCorrection
{
private:
  
  // Dummy constructor
  G4PairingCorrection();
  
  static G4PairingCorrection* theInstance;
   
public:
	
  static G4PairingCorrection* GetInstance();
  
  ~G4PairingCorrection() {};

  G4double GetPairingCorrection(const G4int A, const G4int Z) const
  {
    G4double PCorrection = 0.0;
    const G4int N = A - Z;
    if (theCookPairingCorrections->IsInTableThisN(N) &&
	theCookPairingCorrections->IsInTableThisZ(Z)) 
      PCorrection = theCookPairingCorrections->GetParingCorrection(A,Z);
    else if (theCameronGilbertPairingCorrections->IsInTableThisN(N) &&
	     theCameronGilbertPairingCorrections->IsInTableThisZ(Z))
      PCorrection = theCameronGilbertPairingCorrections->GetPairingCorrection(A,Z);
    else {
      const G4double PairingConstant = 12.0*MeV;
      G4double Pair = (1.0 - static_cast<G4double>(Z) + 2.0*(Z/2)) + (1.0 - static_cast<G4double>(N) + 2.0*(N/2));
      PCorrection = Pair*PairingConstant/sqrt(static_cast<G4double>(A));
    }
    return std::max(PCorrection,0.0);
  }


  G4double GetFissionPairingCorrection(const G4int A, const G4int Z) const 
  {
    const G4double PairingConstant = 14.0*MeV;
    const G4int N = A - Z;
    G4double Pair = (1.0 - static_cast<G4double>(Z) + 2.0*(Z/2)) + (1.0 - static_cast<G4double>(N) + 2.0*(N/2));
    G4double PCorrection = Pair*PairingConstant/sqrt(static_cast<G4double>(A));
    return PCorrection;
  }

private:

  
  G4CookPairingCorrections* theCookPairingCorrections;
//  G4CameronTruranHilfPairingCorrections* theCameronTruranHilfPairingCorrections;
  G4CameronGilbertPairingCorrections* theCameronGilbertPairingCorrections;

};
#endif
