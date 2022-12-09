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
/*
 * G4O3.cc
 *
 *  Created on: Jul 23, 2019
 *      Author: W. G. Shin
 *              J. Ramos-Mendez and B. Faddegon
*/

#include "G4O3.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"

// ######################################################################
// ###                         Ozone                               ###
// ######################################################################
G4O3* G4O3::theInstance = 0;

G4O3* G4O3::Definition()
{
  if (theInstance != 0) return theInstance;
  const G4String name = "O_3";
  // search in particle table]
  G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* anInstance = pTable->FindParticle(name);
  if (anInstance == 0)
  {
    const G4String formatedName = "O_{3}";

    // create molecule
    //
    //     G4MoleculeDefinition(const G4String& name,
    //                          G4double mass,
    //                          G4double diffCoeff,
    //                          G4int 	 charge = 0,
    //                          G4int    electronicLevels = 0,
    //                          G4double radius = -1,
    //                          G4int    atomsNumber = -1,
    //                          G4double lifetime = -1,
    //                          G4String aType = "",
    //                          G4FakeParticleID ID = G4FakeParticleID::Create()
    //    );

    G4double mass = 47.99319 * g / Avogadro * c_squared;
    anInstance = new G4MoleculeDefinition(name, mass, 2.0e-9 * (m * m / s), 0, 0,
                                          2.0 * angstrom, 3);

    ((G4MoleculeDefinition*) anInstance)->SetLevelOccupation(0); // Set 2 electrons on 1 single occupancy
    ((G4MoleculeDefinition*) anInstance)->SetFormatedName(formatedName);
  }
  theInstance = static_cast<G4O3*>(anInstance);
  return theInstance;
}
