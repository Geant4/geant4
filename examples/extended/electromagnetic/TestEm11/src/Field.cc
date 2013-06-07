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
/// \file electromagnetic/TestEm11/src/Field.cc
/// \brief Implementation of the Field class
//
// $Id: Field.cc 
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......




#include "Field.hh"
#include "FieldMessenger.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
//#include "G4UniformMagField.hh"

//#include "globals.hh"
#include "G4SystemOfUnits.hh"

Field::Field()
  : G4UniformMagField(G4ThreeVector())
{
  fFieldMessenger = new FieldMessenger(this);
  SetMagFieldValue(0);
  G4cout << "Constructor Field() called" << G4endl;
}

Field::Field(G4double val)
  : G4UniformMagField(G4ThreeVector())
{
  fFieldMessenger = new FieldMessenger(this);
  SetMagFieldValue(val);
  G4cout << "Constructor Field(G4double) called" << G4endl;
}

Field::~Field()
{
  delete fFieldMessenger;
}

//void Field::GetFieldValue(const double*, double* BF) 
//{ BF[0]=fieldVal.x(); BF[1]=fieldVal.y(); BF[2]=fieldVal.z(); }

void Field::SetMagFieldValue(G4double val)
{
  G4cout << "Setting field to [T]: " << val/tesla << G4endl;
  //fFieldVal = val;
  G4FieldManager* fieldMan
    = G4TransportationManager::GetTransportationManager()->GetFieldManager();
  if (val) {
    SetFieldValue(G4ThreeVector(0.,0.,val));
    fieldMan->SetDetectorField(this);
    fieldMan->CreateChordFinder(this);
  }
  else {
    fieldMan->SetDetectorField(NULL);
  }
}

