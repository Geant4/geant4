// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: F03ElectroMagneticField.cc,v 1.1 2001-06-08 11:55:57 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//  
//   User Field class implementation.
//


#include "F03ElectroMagneticField.hh"
#include "F03FieldMessenger.hh"

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

F03ElectroMagneticField::F03ElectroMagneticField()
  :  fStepper(NULL),fChordFinder(NULL),fLocalChordFinder(NULL)
{
  fMagneticField = new G4UniformMagField(
		       G4ThreeVector(3.3*tesla,
                                     0.0,         // 0.5*tesla,
                                     0.0               ));
  fLocalMagneticField = new G4UniformMagField(
		            G4ThreeVector(3.3*tesla,
                                          0.0,         // 0.5*tesla,
                                          0.0               ));

  //  G4FieldManager* fieldMgr = G4TransportationManager::GetTransportationManager()
  //                         ->GetFieldManager();  
  //  fieldMgr->SetDetectorField(fMagneticField);
  //  fieldMgr->CreateChordFinder(fMagneticField);

  fFieldMessenger = new F03FieldMessenger(this) ;  
 
  fEquation = new G4Mag_UsualEqRhs(fMagneticField); 
  fLocalEquation = new G4Mag_UsualEqRhs(fLocalMagneticField); 
 
  fMinStep     = 1.0*mm ; // minimal step of 1 mm is default

  fStepperType = 4 ;      // ClassicalRK4 is default stepper

  fFieldManager = G4TransportationManager::GetTransportationManager()
                                         ->GetFieldManager();

  fLocalFieldManager = G4TransportationManager::GetTransportationManager()
                                         ->GetFieldManager();

  UpdateField();

  //  GetGlobalFieldManager()->CreateChordFinder(this);
}

/////////////////////////////////////////////////////////////////////////////////

F03ElectroMagneticField::F03ElectroMagneticField(G4ThreeVector fieldVector)
{    
  fMagneticField = new G4UniformMagField(fieldVector);
  GetGlobalFieldManager()->CreateChordFinder(this);
}

////////////////////////////////////////////////////////////////////////////////

F03ElectroMagneticField::~F03ElectroMagneticField()
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

void F03ElectroMagneticField::UpdateField()
{
  SetStepper();
  G4cout<<"The minimal step is equal to "<<fMinStep/mm<<" mm"<<G4endl ;

  fFieldManager->SetDetectorField(fMagneticField );
  fLocalFieldManager->SetDetectorField(fLocalMagneticField );

  if(fChordFinder) delete fChordFinder;
  if(fLocalChordFinder) delete fLocalChordFinder;

  fChordFinder = new G4ChordFinder( fMagneticField, fMinStep,fStepper);
  fLocalChordFinder = new G4ChordFinder( fLocalMagneticField, 
                                         fMinStep,fLocalStepper);

  fFieldManager->SetChordFinder( fChordFinder );
  fLocalFieldManager->SetChordFinder( fLocalChordFinder );

  return;
}

/////////////////////////////////////////////////////////////////////////////
//
// Set stepper according to the stepper type
//

void F03ElectroMagneticField::SetStepper()
{
  if(fStepper) delete fStepper;

  switch ( fStepperType ) 
  {
    case 0:  
      fStepper = new G4ExplicitEuler( fEquation ); 
      fLocalStepper = new G4ExplicitEuler( fLocalEquation ); 
      G4cout<<"G4ExplicitEuler is calledS"<<G4endl;     
      break;
    case 1:  
      fStepper = new G4ImplicitEuler( fEquation );      
      fLocalStepper = new G4ImplicitEuler( fLocalEquation );      
      G4cout<<"G4ImplicitEuler is called"<<G4endl;     
      break;
    case 2:  
      fStepper = new G4SimpleRunge( fEquation );        
      fLocalStepper = new G4SimpleRunge( fLocalEquation );        
      G4cout<<"G4SimpleRunge is called"<<G4endl;     
      break;
    case 3:  
      fStepper = new G4SimpleHeum( fEquation );         
      fLocalStepper = new G4SimpleHeum( fLocalEquation );         
      G4cout<<"G4SimpleHeum is called"<<G4endl;     
      break;
    case 4:  
      fStepper = new G4ClassicalRK4( fEquation );       
      fLocalStepper = new G4ClassicalRK4( fLocalEquation );       
      G4cout<<"G4ClassicalRK4 (default) is called"<<G4endl;     
      break;
    case 5:  
      fStepper = new G4HelixExplicitEuler( fEquation ); 
      fLocalStepper = new G4HelixExplicitEuler( fLocalEquation ); 
      G4cout<<"G4HelixExplicitEuler is called"<<G4endl;     
      break;
    case 6:  
      fStepper = new G4HelixImplicitEuler( fEquation ); 
      fLocalStepper = new G4HelixImplicitEuler( fLocalEquation ); 
      G4cout<<"G4HelixImplicitEuler is called"<<G4endl;     
      break;
    case 7:  
      fStepper = new G4HelixSimpleRunge( fEquation );   
      fLocalStepper = new G4HelixSimpleRunge( fLocalEquation );   
      G4cout<<"G4HelixSimpleRunge is called"<<G4endl;     
      break;
    case 8:  
      fStepper = new G4CashKarpRKF45( fEquation );      
      fLocalStepper = new G4CashKarpRKF45( fLocalEquation );      
      G4cout<<"G4CashKarpRKF45 is called"<<G4endl;     
      break;
    case 9:  
      fStepper = new G4RKG3_Stepper( fEquation );       
      fLocalStepper = new G4RKG3_Stepper( fLocalEquation );       
      G4cout<<"G4RKG3_Stepper is called"<<G4endl;     
      break;
    default: fStepper = 0;
  }
  return; 
}

/////////////////////////////////////////////////////////////////////////////
//
// Set the value of the Global Field to fieldValue along Z
//

void F03ElectroMagneticField::SetFieldValue(G4double fieldValue)
{
  if(fMagneticField) delete fMagneticField;
  fMagneticField = new  G4UniformMagField(G4ThreeVector(0,0,fieldValue));
  //  UpdateField();
}

///////////////////////////////////////////////////////////////////////////////
//
// Set the value of the Global Field
//

void F03ElectroMagneticField::SetFieldValue(G4ThreeVector fieldVector)
{
  // Find the Field Manager for the global field
  G4FieldManager* fieldMgr= GetGlobalFieldManager();
    
  if(fieldVector != G4ThreeVector(0.,0.,0.))
  { 
    if(fMagneticField) delete fMagneticField;
    fMagneticField = new  G4UniformMagField(fieldVector);

    // UpdateField();
   
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

G4FieldManager*  F03ElectroMagneticField::GetGlobalFieldManager()
{
  return G4TransportationManager::GetTransportationManager()
	                        ->GetFieldManager();
}
    



