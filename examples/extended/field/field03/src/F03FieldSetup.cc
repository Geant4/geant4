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
// $Id: F03FieldSetup.cc,v 1.2 2003/12/01 17:29:57 japost Exp $
// GEANT4 tag $Name: geant4-06-00 $
//
//  
//   Field Setup class implementation.
//

#include "F03FieldSetup.hh"
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

F03FieldSetup::F03FieldSetup()
  :  fChordFinder(0), fLocalChordFinder(0), fStepper(0)
{
  fMagneticField = new G4UniformMagField(
		       G4ThreeVector(3.3*tesla,
                                     0.0,              // 0.5*tesla,
                                     0.0       ));
  fLocalMagneticField = new G4UniformMagField(
		            G4ThreeVector(3.3*tesla,
                                          0.0,         // 0.5*tesla,
                                          0.0  ));

  fFieldMessenger = new F03FieldMessenger(this) ;  
 
  fEquation = new G4Mag_UsualEqRhs(fMagneticField); 
  fLocalEquation = new G4Mag_UsualEqRhs(fLocalMagneticField); 
 
  fMinStep     = 0.25*mm ; // minimal step of 1 mm is default
  fStepperType = 4 ;      // ClassicalRK4 is default stepper

  fFieldManager = fLocalFieldManager = GetGlobalFieldManager();

  UpdateField();
}

/////////////////////////////////////////////////////////////////////////////////

F03FieldSetup::F03FieldSetup(G4ThreeVector fieldVector)
{    
  fMagneticField = new G4UniformMagField(fieldVector);
  GetGlobalFieldManager()->CreateChordFinder(fMagneticField);
}

////////////////////////////////////////////////////////////////////////////////

F03FieldSetup::~F03FieldSetup()
{
  if(fMagneticField) delete fMagneticField;
  if(fChordFinder)   delete fChordFinder;
  if(fStepper)       delete fStepper;
}

/////////////////////////////////////////////////////////////////////////////
//
// Update field
//

void F03FieldSetup::UpdateField()
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
}

/////////////////////////////////////////////////////////////////////////////
//
// Set stepper according to the stepper type
//

void F03FieldSetup::SetStepper()
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
}

/////////////////////////////////////////////////////////////////////////////
//
// Set the value of the Global Field to fieldValue along Z
//

void F03FieldSetup::SetFieldValue(G4double fieldValue)
{
  if(fMagneticField) delete fMagneticField;
  fMagneticField = new  G4UniformMagField(G4ThreeVector(0,0,fieldValue));
  //  UpdateField();
}

///////////////////////////////////////////////////////////////////////////////
//
// Set the value of the Global Field
//

void F03FieldSetup::SetFieldValue(G4ThreeVector fieldVector)
{
  // Find the Field Manager for the global field
  G4FieldManager* fieldMgr= GetGlobalFieldManager();
    
  if(fMagneticField) delete fMagneticField;

  if(fieldVector != G4ThreeVector(0.,0.,0.))
  { 
    fMagneticField = new  G4UniformMagField(fieldVector);
  }
  else 
  {
    // If the new field's value is Zero, then it is best to
    //  insure that it is not used for propagation.
    fMagneticField = 0; 
  }
  fieldMgr->SetDetectorField(fMagneticField);

  // UpdateField();
}

////////////////////////////////////////////////////////////////////////////////
//
//  Utility method

G4FieldManager*  F03FieldSetup::GetGlobalFieldManager()
{
  return G4TransportationManager::GetTransportationManager()
	                        ->GetFieldManager();
}

