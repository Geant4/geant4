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
// $Id: F02ElectroMagneticField.cc,v 1.2 2001-11-13 17:22:39 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//  
//   User Field class implementation.
//


#include "F02ElectroMagneticField.hh"
#include "F02FieldMessenger.hh"

#include "G4UniformElectricField.hh"
#include "G4UniformMagField.hh"
#include "G4MagneticField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4EquationOfMotion.hh"
#include "G4EqMagElectricField.hh"
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
  :  fStepper(NULL),fChordFinder(NULL),G4UniformElectricField(
                 G4ThreeVector(0.0,1000.0*kilovolt/cm,1000.0*kilovolt/cm))
{
  fEMfield = new G4UniformElectricField(
                 G4ThreeVector(
                                0.0*kilovolt/cm, 
                                0.0*kilovolt/cm,
                                1000.0*kilovolt/cm ) );

  fFieldMessenger = new F02FieldMessenger(this) ;  
 
  fEquation = new G4EqMagElectricField(fEMfield); 
 
  fMinStep     = 1.0*mm ; // minimal step of 1 mm is default

  fStepperType = 4 ;      // ClassicalRK4 is default stepper

  fFieldManager = G4TransportationManager::GetTransportationManager()
                                         ->GetFieldManager();

  UpdateField();

  //  GetGlobalFieldManager()->CreateChordFinder(this);
}

/////////////////////////////////////////////////////////////////////////////////

F02ElectroMagneticField::F02ElectroMagneticField(G4ThreeVector pFV):
G4UniformElectricField(pFV)
{    
  fEMfield = new G4UniformElectricField(pFV);
  GetGlobalFieldManager()->CreateChordFinder(this);
}

////////////////////////////////////////////////////////////////////////////////

F02ElectroMagneticField::~F02ElectroMagneticField()
{
  // GetGlobalFieldManager()->SetDetectorField(0);

  if(fEMfield) delete fEMfield;
  if(fChordFinder)   delete fChordFinder;
  if(fStepper)       delete fStepper;
}

/////////////////////////////////////////////////////////////////////////////
//
// Update field
//

void F02ElectroMagneticField::UpdateField()
{
  SetStepper();
  G4cout<<"The minimal step is equal to "<<fMinStep/mm<<" mm"<<G4endl ;

  fFieldManager->SetDetectorField(fEMfield );

  if(fChordFinder) delete fChordFinder;

  fChordFinder = new G4ChordFinder( fEMfield, fMinStep,fStepper);

  fFieldManager->SetChordFinder( fChordFinder );

  return;
}

/////////////////////////////////////////////////////////////////////////////
//
// Set stepper according to the stepper type
//

void F02ElectroMagneticField::SetStepper()
{
  G4int nvar = 8;

  if(fStepper) delete fStepper;

  switch ( fStepperType ) 
  {
    case 0:  
      fStepper = new G4ExplicitEuler( fEquation, nvar ); 
      G4cout<<"G4ExplicitEuler is calledS"<<G4endl;     
      break;
    case 1:  
      fStepper = new G4ImplicitEuler( fEquation, nvar );      
      G4cout<<"G4ImplicitEuler is called"<<G4endl;     
      break;
    case 2:  
      fStepper = new G4SimpleRunge( fEquation, nvar );        
      G4cout<<"G4SimpleRunge is called"<<G4endl;     
      break;
    case 3:  
      fStepper = new G4SimpleHeum( fEquation, nvar );         
      G4cout<<"G4SimpleHeum is called"<<G4endl;     
      break;
    case 4:  
      fStepper = new G4ClassicalRK4( fEquation, nvar );       
      G4cout<<"G4ClassicalRK4 (default) is called"<<G4endl;     
      break;
    case 5:  
      fStepper = new G4CashKarpRKF45( fEquation, nvar );      
      G4cout<<"G4CashKarpRKF45 is called"<<G4endl;     
      break;
      /* *********************************** 
    case 6:  
      fStepper = new G4RKG3_Stepper( fEquation );       
      G4cout<<"G4RKG3_Stepper is called"<<G4endl;     
      break;
   case 7:  
      fStepper = new G4HelixExplicitEuler( fEquation ); 
      G4cout<<"G4HelixExplicitEuler is called"<<G4endl;     
      break;
    case 8:  
      fStepper = new G4HelixImplicitEuler( fEquation ); 
      G4cout<<"G4HelixImplicitEuler is called"<<G4endl;     
      break;
    case 9:  
      fStepper = new G4HelixSimpleRunge( fEquation );   
      G4cout<<"G4HelixSimpleRunge is called"<<G4endl;     
      break;
      ****************************************** */
    default: fStepper = 0;
  }
  return; 
}

/////////////////////////////////////////////////////////////////////////////
//
// Set the value of the Global Field to fieldValue along Z
//

void F02ElectroMagneticField::SetFieldValue(G4double fieldValue)
{
  if(fEMfield) delete fEMfield;
  fEMfield = new  G4UniformElectricField(G4ThreeVector(0,0,fieldValue));
  //  UpdateField();
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
    if(fEMfield) delete fEMfield;
    fEMfield = new  G4UniformElectricField(fieldVector);

    // UpdateField();
   
    fieldMgr->SetDetectorField(this);
  }
  else 
  {
    // If the new field's value is Zero, then it is best to
    //  insure that it is not used for propagation.

    G4MagneticField* fEMfield = NULL;
    fieldMgr->SetDetectorField(fEMfield);
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
    









