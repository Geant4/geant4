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
// $Id: G4MonopoleFieldSetup.cc,v 1.1 2010-06-04 19:03:36 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

// #include "G4SIunits.hh"

G4MonopoleFieldSetup* G4MonopoleFieldSetup::fMonopoleFieldSetup=0;

//////////////////////////////////////////////////////////////////////////
//
//  Constructor
G4MonopoleFieldSetup::G4MonopoleFieldSetup()
{
  if (!fMonopoleFieldSetup)
  {
    fChordFinder = 0;
    fStepper = 0;
    fMagneticField = new G4UniformMagField( G4ThreeVector( 0.0, 0.0, 3.0*tesla) ); 
    InitialiseAll();
  }
  else
  {
    G4cerr << "Only ONE instance of G4MonopoleFieldSetup is allowed!"
           << G4endl;
    G4Exception("G4MonopoleFieldSetup::G4MonopoleFieldSetup()",
                "InvalidSetup", FatalException,
                "Only ONE instance of MagneticFieldSetup is allowed!");
  }
}

//////////////////////////////////////////////////////////////////////////
// Retrieve the static instance of the singleton
//
G4MonopoleFieldSetup* G4MonopoleFieldSetup::GetMonopoleFieldSetup()
{
   static G4MonopoleFieldSetup theInstance;
   if (!fMonopoleFieldSetup)
   {
     fMonopoleFieldSetup = &theInstance;
   }
   
   return fMonopoleFieldSetup;
}

void
G4MonopoleFieldSetup::InitialiseAll()
{
  fEquation = new G4Mag_UsualEqRhs(fMagneticField); 
  fMonopoleEquation = new G4MonopoleEquation(fMagneticField);
 
  fMinStep     = 0.01*mm ; // minimal step of 1 mm is default

  fFieldManager = G4TransportationManager::GetTransportationManager()
                                         ->GetFieldManager();
                                         
  fMonopoleStepper = new G4ClassicalRK4( fMonopoleEquation );       
  fStepper = new G4ClassicalRK4( fEquation );                                         
                                         
  SetStepperAndChordFinder(0);
}

////////////////////////////////////////////////////////////////////////////////

G4MonopoleFieldSetup::~G4MonopoleFieldSetup()
{
  if(fMagneticField) delete fMagneticField;
  if(fChordFinder)   delete fChordFinder;
  if(fStepper)       delete fStepper;
  if(fMonopoleStepper)  delete fMonopoleStepper;
}

/////////////////////////////////////////////////////////////////////////////
//
// Update field
//
void G4MonopoleFieldSetup::SetStepperAndChordFinder(G4int val)
{
  switch (val)
  {
    case 0:
      fStepper = new G4ClassicalRK4( fEquation );       
//      G4cout<<"G4ClassicalRK4 (default) is called"<<G4endl;     
      break;
    case 1: 
      fStepper = new G4ClassicalRK4( fMonopoleEquation );       
//      G4cout<<"G4ClassicalRK4 (default) is called"<<G4endl;     
      break;
  }   
  
//  G4cout<<"The minimal step is equal to "<<fMinStep/mm<<" mm"<<G4endl ;

  fFieldManager->SetDetectorField(fMagneticField );

  if(fChordFinder) delete fChordFinder;

  fChordFinder = new G4ChordFinder( fMagneticField, fMinStep, fStepper);

  fFieldManager->SetChordFinder( fChordFinder );

  return;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Utility method
G4FieldManager*  G4MonopoleFieldSetup::GetGlobalFieldManager()
{
  return G4TransportationManager::GetTransportationManager()
                          ->GetFieldManager();
}
