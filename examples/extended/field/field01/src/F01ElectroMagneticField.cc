// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: F01ElectroMagneticField.cc,v 1.3 2001-03-29 10:51:00 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//  
//   User Field class implementation.
//


#include "F01ElectroMagneticField.hh"
#include "F01FieldMessenger.hh"

#include "G4UniformMagField.hh"
#include "G4MagneticField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4ChordFinder.hh"

#include "G4ExplicitEuler.hh"
#include "G4ImplicitEuler.hh"
#include "G4SimpleRunge.hh"
#include "G4SimpleHeum.hh"
#include "G4ClassicalRK4.hh"
#include "G4HelixExplicitEuler.hh"
#include "G4HelixImplicitEuler.hh"
#include "G4HelixSimpleRunge.hh"
#include "G4CashKarpRKF45.hh"
#include "G4RKG3_Stepper.hh"

//////////////////////////////////////////////////////////////////////////
//
//  Constructors:

F01ElectroMagneticField::F01ElectroMagneticField()
{
  fMagneticField = new G4UniformMagField(G4ThreeVector(0.0,0.2*tesla,0.0));

  //  G4FieldManager* fieldMgr = G4TransportationManager::GetTransportationManager()
  //                         ->GetFieldManager();  
  //  fieldMgr->SetDetectorField(fMagneticField);
  //  fieldMgr->CreateChordFinder(fMagneticField);

  fFieldMessenger = new F01FieldMessenger(this) ;  
 
  fEquation = new G4Mag_UsualEqRhs(fMagneticField); 
 
  fMinStep     = 1.0*mm ;

  fStepperType = 4 ;

  fFieldManager = G4TransportationManager::GetTransportationManager()
                                         ->GetFieldManager();

  UpdateField();

  //  GetGlobalFieldManager()->CreateChordFinder(this);
}

/////////////////////////////////////////////////////////////////////////////////

F01ElectroMagneticField::F01ElectroMagneticField(G4ThreeVector fieldVector)
{    
  fMagneticField = new G4UniformMagField(fieldVector);
  GetGlobalFieldManager()->CreateChordFinder(this);
}

////////////////////////////////////////////////////////////////////////////////

F01ElectroMagneticField::~F01ElectroMagneticField()
{
  // GetGlobalFieldManager()->SetDetectorField(0);

  if(fMagneticField) delete fMagneticField;
  if(fChordFinder)   delete fChordFinder;
  if(fStepper)       delete fStepper;
}

/////////////////////////////////////////////////////////////////////////////
//
// Update field
//

void F01ElectroMagneticField::UpdateField()
{
  SetStepper();

  fFieldManager->SetDetectorField(fMagneticField );

  if(fChordFinder) delete fChordFinder;
  fChordFinder = new G4ChordFinder( fMagneticField, fMinStep,fStepper);

  fFieldManager->SetChordFinder( fChordFinder );

  return;
}

/////////////////////////////////////////////////////////////////////////////
//
// Set stepper according to the stepper type
//

void F01ElectroMagneticField::SetStepper()
{
  if(fStepper) delete fStepper;

  switch ( fStepperType ) 
  {
    case 0:  fStepper = new G4ExplicitEuler( fEquation );      break;
    case 1:  fStepper = new G4ImplicitEuler( fEquation );      break;
    case 2:  fStepper = new G4SimpleRunge( fEquation );        break;
    case 3:  fStepper = new G4SimpleHeum( fEquation );         break;
    case 4:  fStepper = new G4ClassicalRK4( fEquation );       break;
    case 5:  fStepper = new G4HelixExplicitEuler( fEquation ); break;
    case 6:  fStepper = new G4HelixImplicitEuler( fEquation ); break;
    case 7:  fStepper = new G4HelixSimpleRunge( fEquation );   break;
    case 8:  fStepper = new G4CashKarpRKF45( fEquation );      break;
    case 9:  fStepper = new G4RKG3_Stepper( fEquation );       break;
    default: fStepper = 0;
  }
  return; 
}

/////////////////////////////////////////////////////////////////////////////
//
// Set the value of the Global Field to fieldValue along Z
//

void F01ElectroMagneticField::SetFieldValue(G4double fieldValue)
{
  if(fMagneticField) delete fMagneticField;
  fMagneticField = new  G4UniformMagField(G4ThreeVector(0,0,fieldValue));
  UpdateField();
}

///////////////////////////////////////////////////////////////////////////////
//
// Set the value of the Global Field
//

void F01ElectroMagneticField::SetFieldValue(G4ThreeVector fieldVector)
{
  // Find the Field Manager for the global field
  G4FieldManager* fieldMgr= GetGlobalFieldManager();
    
  if(fieldVector != G4ThreeVector(0.,0.,0.))
  { 
    if(fMagneticField) delete fMagneticField;
    fMagneticField = new  G4UniformMagField(fieldVector);

    UpdateField();
   
    fieldMgr->SetDetectorField(this);
  }
  else 
  {
    // If the new field's value is Zero, then it is best to
    //  insure that it is not used for propagation.

    G4MagneticField* fMagneticField = NULL;
    fieldMgr->SetDetectorField(fMagneticField);
  }
}

////////////////////////////////////////////////////////////////////////////////
//
//  Utility method

G4FieldManager*  F01ElectroMagneticField::GetGlobalFieldManager()
{
  return G4TransportationManager::GetTransportationManager()
	                        ->GetFieldManager();
}
    



