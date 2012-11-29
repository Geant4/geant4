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
// $Id$
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara

#include "G4PairingCorrection.hh"
#include "G4SystemOfUnits.hh"

G4PairingCorrection* G4PairingCorrection::theInstance = 0;

G4PairingCorrection::G4PairingCorrection()
{
  theCookPairingCorrections =  G4CookPairingCorrections::GetInstance();
  theCameronGilbertPairingCorrections = G4CameronGilbertPairingCorrections::GetInstance();
//  theCameronTruranHilfPairingCorrections = G4CameronTruranHilfPairingCorrections::GetInstance();
}

G4PairingCorrection::~G4PairingCorrection()
{;}

G4PairingCorrection* G4PairingCorrection::GetInstance()
{
  if (!theInstance)  { 
    static G4PairingCorrection theCorrections;
    theInstance = &theCorrections; 
  }
  return theInstance;
}   

G4double G4PairingCorrection::GetPairingCorrection(G4int A, G4int Z) const
{
  G4double PCorrection = 0.0;
  G4int N = A - Z;
  if (theCookPairingCorrections->IsInTableThisN(N) &&
      theCookPairingCorrections->IsInTableThisZ(Z)) 
    PCorrection = theCookPairingCorrections->GetParingCorrection(A,Z);
  else if (theCameronGilbertPairingCorrections->IsInTableThisN(N) &&
	   theCameronGilbertPairingCorrections->IsInTableThisZ(Z))
    PCorrection = theCameronGilbertPairingCorrections->GetPairingCorrection(A,Z);
  else {
    const G4double PairingConstant = 12.0*MeV;
    G4double Pair = (1 - Z + 2*(Z/2)) + (1 - N + 2*(N/2));
    PCorrection = Pair*PairingConstant/std::sqrt(static_cast<G4double>(A));
  }
  return std::max(PCorrection,0.0);
}


G4double G4PairingCorrection::GetFissionPairingCorrection(G4int A, G4int Z) const 
{
  const G4double PairingConstant = 14.0*MeV;
  G4int N = A - Z;
  G4double Pair = (1 - Z + 2*(Z/2)) + (1 - N + 2*(N/2));
  G4double PCorrection = Pair*PairingConstant/std::sqrt(static_cast<G4double>(A));
  return PCorrection;
}
