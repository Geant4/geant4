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
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4PhotoElectricAngularGenerator462
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

#include "G4PhotoElectricAngularGenerator462.hh"
#include "Randomize.hh"

//    

G4PhotoElectricAngularGenerator462::G4PhotoElectricAngularGenerator462(const G4String& name):G4VPhotoElectricAngularDistribution(name)
{;}

//    

G4PhotoElectricAngularGenerator462::~G4PhotoElectricAngularGenerator462() 
{;}

//

G4ThreeVector G4PhotoElectricAngularGenerator462::GetPhotoElectronDirection(G4ThreeVector direction, G4double)
{
  return direction;
}

//

void G4PhotoElectricAngularGenerator462::PrintGeneratorInformation() const
{
  G4cout << "\n" << G4endl;
  G4cout << "Simple Photoelectric Angular Generator" << G4endl;
  G4cout << "Photoelectron is emmited with the same direction " << G4endl;
  G4cout << "than the incident photon (see Physics Reference Manual) \n" << G4endl;
} 
