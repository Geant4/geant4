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
// $Id: F03FieldSetup.cc,v 1.5 2007-04-28 01:31:12 gum Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

  fFieldManager = GetGlobalFieldManager();
  fLocalFieldManager = new G4FieldManager();

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

void F03FieldSetup::SetFieldValue(G4double fieldStrength)
{
  G4ThreeVector fieldSetVec(0.0, 0.0, fieldStrength);
  this->SetFieldValue( fieldSetVec ); 
  //    *************

}

///////////////////////////////////////////////////////////////////////////////
//
// Set the value of the Global Field
//

void F03FieldSetup::SetFieldValue(G4ThreeVector fieldVector)
{
  if(fMagneticField) delete fMagneticField;

  if(fieldVector != G4ThreeVector(0.,0.,0.))
  { 
    fMagneticField = new  G4UniformMagField(fieldVector);
  }
  else 
  {
    // If the new field's value is Zero, then 
    //  setting the pointer to zero ensures 
    //  that it is not used for propagation.
    fMagneticField = 0; 
  }

  // Either  
  //   - UpdateField() to reset all (ChordFinder, Equation);
     // UpdateField();
  //   or simply update the field manager & equation of motion 
  //      with pointer to new field
  GetGlobalFieldManager()->SetDetectorField(fMagneticField);
  fEquation->SetFieldObj( fMagneticField ); 

}

////////////////////////////////////////////////////////////////////////////////
//
//  Utility method

G4FieldManager*  F03FieldSetup::GetGlobalFieldManager()
{
  return G4TransportationManager::GetTransportationManager()
	                        ->GetFieldManager();
}


// In place of G4UniformField::GetConstantFieldValue ...
// 
G4ThreeVector F03FieldSetup::GetConstantFieldValue()
{
  static G4double fieldValue[6],  position[4]; 
  position[0] = position[1] = position[2] = position[3] = 0.0; 

  fMagneticField->GetFieldValue( position, fieldValue);
  G4ThreeVector fieldVec(fieldValue[0], fieldValue[1], fieldValue[2]); 

  return fieldVec;
}
