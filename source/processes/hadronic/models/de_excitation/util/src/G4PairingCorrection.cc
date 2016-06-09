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
// $Id: G4PairingCorrection.cc,v 1.4.2.1 2009/03/04 14:56:06 gcosmo Exp $
// GEANT4 tag $Name: geant4-09-02-patch-01 $
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

G4PairingCorrection::~G4PairingCorrection()
{ delete theInstance; }

G4PairingCorrection* G4PairingCorrection::GetInstance()
{
  if (!theInstance) theInstance = new G4PairingCorrection();
  return theInstance;
}   
