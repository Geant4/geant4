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
// $Id: G4PairingCorrection.hh 105799 2017-08-21 07:35:55Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara
//
// Modified:
// 21.03.2013 V.Ivanchenko redesigned and cleaned up

#ifndef G4PairingCorrection_h
#define G4PairingCorrection_h 1

#include "globals.hh"
#include "G4CookPairingCorrections.hh"
#include "G4CameronGilbertPairingCorrections.hh"

class G4PairingCorrection
{
public:
	
  explicit G4PairingCorrection();
  
  ~G4PairingCorrection();

  G4double GetPairingCorrection(G4int A, G4int Z) const;

  G4double GetFissionPairingCorrection(G4int A, G4int Z) const;

private:

  G4CookPairingCorrections theCookPairingCorrections;
  G4CameronGilbertPairingCorrections theCameronGilbertPairingCorrections;

};
#endif
