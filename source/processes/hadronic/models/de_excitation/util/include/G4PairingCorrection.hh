//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4PairingCorrection.hh,v 1.6 2009/03/04 11:05:02 gcosmo Exp $
// GEANT4 tag $Name: geant4-09-03 $
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
  
  ~G4PairingCorrection();

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
      PCorrection = Pair*PairingConstant/std::sqrt(static_cast<G4double>(A));
    }
    return std::max(PCorrection,0.0);
  }


  G4double GetFissionPairingCorrection(const G4int A, const G4int Z) const 
  {
    const G4double PairingConstant = 14.0*MeV;
    const G4int N = A - Z;
    G4double Pair = (1.0 - static_cast<G4double>(Z) + 2.0*(Z/2)) + (1.0 - static_cast<G4double>(N) + 2.0*(N/2));
    G4double PCorrection = Pair*PairingConstant/std::sqrt(static_cast<G4double>(A));
    return PCorrection;
  }

private:

  
  G4CookPairingCorrections* theCookPairingCorrections;
//  G4CameronTruranHilfPairingCorrections* theCameronTruranHilfPairingCorrections;
  G4CameronGilbertPairingCorrections* theCameronGilbertPairingCorrections;

};
#endif
