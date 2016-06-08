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
// $Id: G4PairingCorrection.cc,v 1.4 2001/10/05 16:13:43 hpw Exp $
// GEANT4 tag $Name: geant4-04-00 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara


#include "G4PairingCorrection.hh"


// G4double G4PairingCorrection::
// GetPairingCorrection(const G4int anA, const G4int aZ)
// {
//   const G4double PairingConstant = 12.0*MeV;
//   const G4int N = anA - aZ;
//   G4double Pair = (1.0 - G4double(aZ) + 2.0*(aZ/2)) + (1.0 - G4double(N) + 2.0*(N/2));
//   G4double PCorrection = Pair*PairingConstant/sqrt(G4double(anA));
//   return PCorrection;
// }



G4PairingCorrection G4PairingCorrection::theInstance(10.0);

G4PairingCorrection::G4PairingCorrection(G4double dummy)
{
    G4double even_more_dummy = dummy;
    even_more_dummy/=2.;
}
