// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: F02ElectroMagneticField.cc,v 1.1 2001-03-27 16:26:56 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//  
//   User Field class implementation.
//


#include "F02ElectroMagneticField.hh"

#include "G4UniformMagField.hh"
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

F02ElectroMagneticField::F02ElectroMagneticField()
  : G4UniformMagField(G4ThreeVector())
{

  G4UniformMagField* magField = new G4UniformMagField(
                                    G4ThreeVector(0.0,0.2*tesla,0.0));

  //  G4FieldManager* fieldMgr = G4TransportationManager::GetTransportationManager()
  //                         ->GetFieldManager();  
  //  fieldMgr->SetDetectorField(magField);
  //  fieldMgr->CreateChordFinder(magField);

  G4FieldManager*         pFieldMgr;
  G4ChordFinder*          pChordFinder;
  G4Mag_UsualEqRhs*       fEquation = new G4Mag_UsualEqRhs(magField); 
  G4MagIntegratorStepper* pStepper;

  fStepperType = 4 ;

  switch ( fStepperType ) 
  {
    case 0:  pStepper = new G4ExplicitEuler( fEquation );      break;
    case 1:  pStepper = new G4ImplicitEuler( fEquation );      break;
    case 2:  pStepper = new G4SimpleRunge( fEquation );        break;
    case 3:  pStepper = new G4SimpleHeum( fEquation );         break;
    case 4:  pStepper = new G4ClassicalRK4( fEquation );       break;
    case 5:  pStepper = new G4HelixExplicitEuler( fEquation ); break;
    case 6:  pStepper = new G4HelixImplicitEuler( fEquation ); break;
    case 7:  pStepper = new G4HelixSimpleRunge( fEquation );   break;
      // case 8:  pStepper = new G4CashKarpRKF45( fEquation );      break;
    case 9:  pStepper = new G4RKG3_Stepper( fEquation );       break;
    default: pStepper = 0;
  }
  pFieldMgr = G4TransportationManager::GetTransportationManager()->GetFieldManager();

  pFieldMgr->SetDetectorField(magField );

  G4double minStep = 1.0*mm ;

  pChordFinder = new G4ChordFinder( magField, minStep,pStepper);

  pFieldMgr->SetChordFinder( pChordFinder );



  //  GetGlobalFieldManager()->CreateChordFinder(this);
}

/////////////////////////////////////////////////////////////////////////////////

F02ElectroMagneticField::F02ElectroMagneticField(G4ThreeVector fieldVector)
  : G4UniformMagField(fieldVector)
{    
  GetGlobalFieldManager()->CreateChordFinder(this);
}

////////////////////////////////////////////////////////////////////////////////

F02ElectroMagneticField::~F02ElectroMagneticField()
{
  // GetGlobalFieldManager()->SetDetectorField(0);
}

/////////////////////////////////////////////////////////////////////////////
//
// Set the value of the Global Field to fieldValue along Z
//

void F02ElectroMagneticField::SetFieldValue(G4double fieldValue)
{
   G4UniformMagField::SetFieldValue(G4ThreeVector(0,0,fieldValue));
}

///////////////////////////////////////////////////////////////////////////////
//
// Set the value of the Global Field
//

void F02ElectroMagneticField::SetFieldValue(G4ThreeVector fieldVector)
{
  // Find the Field Manager for the global field
  G4FieldManager* fieldMgr= GetGlobalFieldManager();
    
  if(fieldVector != G4ThreeVector(0.,0.,0.))
  { 
    G4UniformMagField::SetFieldValue(fieldVector);
    fieldMgr->SetDetectorField(this);
  }
  else 
  {
    // If the new field's value is Zero, then it is best to
    //  insure that it is not used for propagation.

    G4MagneticField* magField = NULL;
    fieldMgr->SetDetectorField(magField);
  }
}

////////////////////////////////////////////////////////////////////////////////
//
//  Utility method

G4FieldManager*  F02ElectroMagneticField::GetGlobalFieldManager()
{
  return G4TransportationManager::GetTransportationManager()
	     ->GetFieldManager();
}
    



