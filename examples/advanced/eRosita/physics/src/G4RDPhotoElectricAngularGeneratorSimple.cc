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
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4RDPhotoElectricAngularGeneratorSimple
//
// Author:        Andreia Trindade (andreia@lip.pt)
// 
// Creation date: 10 May 2004
//
// Modifications: 
// 10 May 2003       A. Trindade    First implementation acording with new design
//
// Class Description: 
//
// Concrete class for PhotoElectric Electron Angular Distribution Generation - ( < 4.6.2 model)
//
// Class Description: End 
//
// -------------------------------------------------------------------
//
//    

#include "G4RDPhotoElectricAngularGeneratorSimple.hh"
#include "Randomize.hh"

//    

G4RDPhotoElectricAngularGeneratorSimple::G4RDPhotoElectricAngularGeneratorSimple(const G4String& name):G4RDVPhotoElectricAngularDistribution(name)
{;}

//    

G4RDPhotoElectricAngularGeneratorSimple::~G4RDPhotoElectricAngularGeneratorSimple() 
{;}

//

G4ThreeVector G4RDPhotoElectricAngularGeneratorSimple::GetPhotoElectronDirection(const G4ThreeVector& direction, const G4double, const G4ThreeVector&, const G4int) const
{
  return direction;
}

//

void G4RDPhotoElectricAngularGeneratorSimple::PrintGeneratorInformation() const
{
  G4cout << "\n" << G4endl;
  G4cout << "Simple Photoelectric Angular Generator" << G4endl;
  G4cout << "Photoelectron is emmited with the same direction " << G4endl;
  G4cout << "than the incident photon (see Physics Reference Manual) \n" << G4endl;
} 
