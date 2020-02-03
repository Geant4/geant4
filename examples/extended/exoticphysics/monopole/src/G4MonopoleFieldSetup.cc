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
/// \file exoticphysics/monopole/src/G4MonopoleFieldSetup.cc
/// \brief Implementation of the G4MonopoleFieldSetup class
//
<<<<<<< HEAD
// $Id: G4MonopoleFieldSetup.cc 68036 2013-03-13 14:13:45Z gcosmo $
=======
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
//
// G4MonopoleFieldSetup is responsible for setting up a magnetic field
// and the ability to use it with two different equation of motions, 
// one for monopoles and another for the rest of the particles. 
// 
// 

// =======================================================================
// Created:  13 May 2010, B. Bozsogi
// =======================================================================

#include "G4MonopoleFieldSetup.hh"
#include "G4MonopoleFieldMessenger.hh"

#include "G4MagneticField.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4MonopoleEquation.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4ChordFinder.hh"

// #include "G4ExplicitEuler.hh"
// #include "G4ImplicitEuler.hh"
// #include "G4SimpleRunge.hh"
// #include "G4SimpleHeum.hh"
#include "G4ClassicalRK4.hh"
// #include "G4HelixExplicitEuler.hh"
// #include "G4HelixImplicitEuler.hh"
// #include "G4HelixSimpleRunge.hh"
// #include "G4CashKarpRKF45.hh"
// #include "G4RKG3_Stepper.hh"

#include "G4SystemOfUnits.hh"

G4MonopoleFieldSetup* G4MonopoleFieldSetup::fMonopoleFieldSetup=0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4MonopoleFieldSetup::G4MonopoleFieldSetup()
 : fFieldManager(0),
   fChordFinder(0),
   fEquation(0),
   fMonopoleEquation(0),
   fMagneticField(0),
   fStepper(0),
   fMonopoleStepper(0),
   fMinStep(0.0),
<<<<<<< HEAD
=======
   fZmagFieldValue(0.2*tesla),
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
   fMonopoleFieldMessenger(0)
{
  fMonopoleFieldMessenger = new G4MonopoleFieldMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4MonopoleFieldSetup* G4MonopoleFieldSetup::GetMonopoleFieldSetup()
{
   if (0 == fMonopoleFieldSetup)
   {
     static G4MonopoleFieldSetup theInstance;
     fMonopoleFieldSetup = &theInstance;
   }
   
   return fMonopoleFieldSetup;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4MonopoleFieldSetup::~G4MonopoleFieldSetup()
{
  delete fMonopoleFieldMessenger;
  if(fMagneticField) delete fMagneticField;
  if(fChordFinder)   delete fChordFinder;
  if(fStepper)       delete fStepper;
  if(fMonopoleStepper)  delete fMonopoleStepper;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

<<<<<<< HEAD
void G4MonopoleFieldSetup::SetMagField(G4double fieldValue)
{
=======
void G4MonopoleFieldSetup::SetZMagFieldValue (G4double val)
{
  //set new magnetic field value and rebuild the field
  fZmagFieldValue = val;

  ConstructMagField();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4MonopoleFieldSetup::ConstructMagField()
{
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
  //apply a global uniform magnetic field along Z axis  
  if (fMagneticField) { delete fMagneticField; }  //delete the existing magn field

  if (fieldValue != 0.)     // create a new one if non nul
    {
<<<<<<< HEAD
      fMagneticField = new G4UniformMagField(G4ThreeVector(0., 0., fieldValue));        
=======
      // G4cout << "Go to create new field ..." << G4endl;
      fMagneticField = new G4UniformMagField(G4ThreeVector(0., 0., 
                                                           fZmagFieldValue));        
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
      InitialiseAll();
    }
   else
    {
      // G4cout << "Set field = 0 ..." << G4endl;
      fMagneticField = 0;
      fFieldManager->SetDetectorField(fMagneticField);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void
G4MonopoleFieldSetup::InitialiseAll()
{
  fFieldManager = G4TransportationManager::GetTransportationManager()
                                          ->GetFieldManager();
                                          
  fEquation = new G4Mag_UsualEqRhs(fMagneticField); 
  fMonopoleEquation = new G4MonopoleEquation(fMagneticField);
 
  fMinStep     = 0.01*mm ; // minimal step of 1 mm is default
                                        
  fMonopoleStepper = new G4ClassicalRK4( fMonopoleEquation, 8 ); // for time information..
  fStepper = new G4ClassicalRK4( fEquation );                                         
                                         
  SetStepperAndChordFinder(0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4MonopoleFieldSetup::SetStepperAndChordFinder(G4int val)
{
  if (fMagneticField)
  {
    fFieldManager->SetDetectorField(fMagneticField );

    if(fChordFinder) delete fChordFinder;

    switch (val)
    {
      case 0:
        fChordFinder = new G4ChordFinder( fMagneticField, fMinStep, fStepper);      
        break;
      case 1: 
        fChordFinder = new G4ChordFinder( fMagneticField, fMinStep, fMonopoleStepper);
        break;
    }   
  
    fFieldManager->SetChordFinder( fChordFinder );
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4FieldManager*  G4MonopoleFieldSetup::GetGlobalFieldManager()
{
  return G4TransportationManager::GetTransportationManager()
                          ->GetFieldManager();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
