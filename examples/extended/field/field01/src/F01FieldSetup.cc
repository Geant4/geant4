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
// $Id: F01FieldSetup.cc,v 1.7 2006-06-29 17:16:36 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//   User Field setup class implementation.
//

#include "F01FieldSetup.hh"
#include "F01FieldMessenger.hh"

#include "G4MagneticField.hh"
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

// #include "G4SIunits.hh"

//////////////////////////////////////////////////////////////////////////
//
//  Constructors:

F01FieldSetup::F01FieldSetup(G4ThreeVector fieldVector)
  : fChordFinder(0), fStepper(0)
{
  fMagneticField = new G4UniformMagField( fieldVector ); 
    // G4ThreeVector(3.3*tesla, 0.0, 0.0 ));
  G4cout << " F01FieldSetup: magnetic field set to Uniform( "
	 << fieldVector << " ) " << G4endl;
  InitialiseAll();
}

F01FieldSetup::F01FieldSetup()
  : fChordFinder(0), fStepper(0)
{
  fMagneticField = new G4UniformMagField( G4ThreeVector(0.0, 0.0, 0.0 ) );
  G4cout << " F01FieldSetup: magnetic field set to Uniform( 0.0, 0, 0 ) " << G4endl;
  InitialiseAll();
}

void
F01FieldSetup::InitialiseAll()
{
  fFieldMessenger = new F01FieldMessenger(this) ;  
 
  fEquation = new G4Mag_UsualEqRhs(fMagneticField); 
 
  fMinStep     = 1.0*mm ; // minimal step of 1 mm is default

  fStepperType = 4 ;      // ClassicalRK4 is default stepper

  fFieldManager = G4TransportationManager::GetTransportationManager()
                                         ->GetFieldManager();
  CreateStepperAndChordFinder();
}

////////////////////////////////////////////////////////////////////////////////

F01FieldSetup::~F01FieldSetup()
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

void F01FieldSetup::CreateStepperAndChordFinder()
{
  SetStepper();
  G4cout<<"The minimal step is equal to "<<fMinStep/mm<<" mm"<<G4endl ;

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

void F01FieldSetup::SetStepper()
{
  if(fStepper) delete fStepper;

  switch ( fStepperType ) 
  {
    case 0:  
      fStepper = new G4ExplicitEuler( fEquation ); 
      G4cout<<"G4ExplicitEuler is calledS"<<G4endl;     
      break;
    case 1:  
      fStepper = new G4ImplicitEuler( fEquation );      
      G4cout<<"G4ImplicitEuler is called"<<G4endl;     
      break;
    case 2:  
      fStepper = new G4SimpleRunge( fEquation );        
      G4cout<<"G4SimpleRunge is called"<<G4endl;     
      break;
    case 3:  
      fStepper = new G4SimpleHeum( fEquation );         
      G4cout<<"G4SimpleHeum is called"<<G4endl;     
      break;
    case 4:  
      fStepper = new G4ClassicalRK4( fEquation );       
      G4cout<<"G4ClassicalRK4 (default) is called"<<G4endl;     
      break;
    case 5:  
      fStepper = new G4HelixExplicitEuler( fEquation ); 
      G4cout<<"G4HelixExplicitEuler is called"<<G4endl;     
      break;
    case 6:  
      fStepper = new G4HelixImplicitEuler( fEquation ); 
      G4cout<<"G4HelixImplicitEuler is called"<<G4endl;     
      break;
    case 7:  
      fStepper = new G4HelixSimpleRunge( fEquation );   
      G4cout<<"G4HelixSimpleRunge is called"<<G4endl;     
      break;
    case 8:  
      fStepper = new G4CashKarpRKF45( fEquation );      
      G4cout<<"G4CashKarpRKF45 is called"<<G4endl;     
      break;
    case 9:  
      fStepper = new G4RKG3_Stepper( fEquation );       
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

void F01FieldSetup::SetFieldValue(G4double fieldStrength)
{
#ifdef G4VERBOSE
  G4cout << "Setting Field strength to " 
	 << fieldStrength / gauss  << " Gauss."; // << G4endl;
#endif

  G4ThreeVector fieldSetVec(0.0, 0.0, fieldStrength);
  this->SetFieldValue( fieldSetVec ); 
  //    *************

#ifdef G4VERBOSE
  G4double fieldValue[6],  position[4]; 
  position[0] = position[1] = position[2] = position[3] = 0.0; 
  if( fieldStrength != 0.0 ){
    fMagneticField->GetFieldValue( position, fieldValue);
    G4ThreeVector fieldVec(fieldValue[0], fieldValue[1], fieldValue[2]); 
    // G4cout << " fMagneticField is now " << fMagneticField 
    G4cout << " Magnetic field vector is " 
	   << fieldVec / gauss << " G " << G4endl; 
  }else{
    if( fMagneticField == 0 )
      G4cout << " Magnetic field pointer is null." << G4endl;
    else
      G4Exception("F01FieldSetup::SetFieldValue(double)",
		  "IncorrectForZeroField",
		  FatalException,
		  "fMagneticField ptr should be set to 0 for no field."); 
  }
#endif
}

///////////////////////////////////////////////////////////////////////////////
//
// Set the value of the Global Field
//

void F01FieldSetup::SetFieldValue(G4ThreeVector fieldVector)
{
  if(fMagneticField) delete fMagneticField;
    
  if(fieldVector != G4ThreeVector(0.,0.,0.))
  { 
    fMagneticField = new  G4UniformMagField(fieldVector);
    // CreateStepperAndChordFinder();
  }
  else 
  {
    // If the new field's value is Zero, signal it as below
    //   so that it is not used for propagation.
    fMagneticField = 0;
  }

  //  Set this as the field of the global Field Manager
  GetGlobalFieldManager()->SetDetectorField(fMagneticField);

  // Now notify equation of new field
  fEquation->SetFieldObj( fMagneticField ); 

}

////////////////////////////////////////////////////////////////////////////////
//
//  Utility method

G4FieldManager*  F01FieldSetup::GetGlobalFieldManager()
{
  return G4TransportationManager::GetTransportationManager()
	                        ->GetFieldManager();
}
