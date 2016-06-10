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
// $Id: G4Li5FermiFragment.cc 67983 2013-03-13 10:42:03Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1998)

#include "G4Li5FermiFragment.hh"
#include "G4NucleiProperties.hh"

G4Li5FermiFragment::G4Li5FermiFragment(G4int anA, G4int aZ, G4int Pol, G4double ExE)
  : G4UnstableFermiFragment(anA,aZ,Pol,ExE)
{
  // Li5 ----> alpha + proton
  Masses.push_back(G4NucleiProperties::GetNuclearMass(4,2));
  Masses.push_back(CLHEP::proton_mass_c2);

  AtomNum.push_back(4);
  AtomNum.push_back(1);
  
  Charges.push_back(2);
  Charges.push_back(1);
}

G4Li5FermiFragment::~G4Li5FermiFragment()
{}
