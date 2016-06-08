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
// $Id: G4PairingCorrection.cc,v 1.5 2001/11/08 10:11:18 vlara Exp $
// GEANT4 tag $Name: geant4-04-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara


#include "G4PairingCorrection.hh"


G4PairingCorrection* G4PairingCorrection::theInstance = 0;

G4PairingCorrection::G4PairingCorrection()
{
  theCookPairingCorrections =  G4CookPairingCorrections::GetInstance();
  theCameronGilbertPairingCorrections = G4CameronGilbertPairingCorrections::GetInstance();
//  theCameronTruranHilfPairingCorrections = G4CameronTruranHilfPairingCorrections::GetInstance();
}

G4PairingCorrection* G4PairingCorrection::GetInstance()
{
  if (!theInstance) theInstance = new G4PairingCorrection();
  return theInstance;
}   
