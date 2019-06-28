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
/// \file field/field01/src/F01FieldSetup.cc
/// \brief Implementation of the F01FieldSetup class
//
//
//
//   User Field setup class implementation.
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

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
#include "G4ConstRK4.hh"
#include "G4NystromRK4.hh"
#include "G4HelixMixedStepper.hh"
#include "G4ExactHelixStepper.hh"

// Newest steppers - from Release 10.3-beta (June 2013)
#include "G4BogackiShampine23.hh"
#include "G4BogackiShampine45.hh"
#include "G4DormandPrince745.hh"
#include "G4DormandPrinceRK56.hh"
#include "G4DormandPrinceRK78.hh"
#include "G4TsitourasRK45.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

enum EStepperNumber { kDormandPrince45 = 17, kBogackiShampine45= 45, kClassicalRK4 = 4,
                      kNystromRK4 = 13 /*soon 40*/,
                      kDormandPrince56 = 56, kBogackiShampine23= 23, kCashKarp = 8,
                      kDormandPrince78 = 78, kTsitouras45 = 145
} ;

//  Constructors:

F01FieldSetup::F01FieldSetup(G4ThreeVector fieldVector,
                             G4int         stepperNum,
                             G4bool        useFSALstepper )
 : fMagneticField(new G4UniformMagField(fieldVector)),
   fUseFSALstepper(useFSALstepper),
   fStepperType(0),
   fMinStep(0.)
{
  G4cout << " F01FieldSetup: magnetic field set to Uniform( "
         << fieldVector << " ) " << G4endl;

  if( stepperNum == -1000 )
  {
     fUseFSALstepper = useFSALstepper;
     if( !useFSALstepper )
        fStepperType=   17;   // Use Dormand Prince (7) 4/5 as default stepper
     else
        fStepperType = 101;
  }
  else
  {
     fUseFSALstepper = ( stepperNum > 0 );
     if( stepperNum > 0 )
        fStepperType =   stepperNum;
     else
        fStepperType = - stepperNum;        
     
  }
  
  InitialiseAll();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F01FieldSetup::F01FieldSetup()
 : fMagneticField(new G4UniformMagField(G4ThreeVector())),
   fUseFSALstepper(false),
   fStepperType(17),   // Use Dormand Prince (7) 4/5 as default stepper
   fMinStep(0.)
{
  G4cout << " F01FieldSetup: magnetic field set to Uniform( 0.0, 0, 0 ) "
         << G4endl;
  InitialiseAll();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F01FieldSetup::InitialiseAll()
{
  fFieldMessenger = new F01FieldMessenger(this);
 
  fEquation = new G4Mag_UsualEqRhs(fMagneticField);
 
  fMinStep     = 1.0*mm; // minimal step of 1 mm is default

  fFieldManager = G4TransportationManager::GetTransportationManager()
                    ->GetFieldManager();

  if( fUseFSALstepper ) {
    CreateFSALStepperAndChordFinder();
  }
  else
  {
    CreateStepperAndChordFinder();
  }

  G4cout  << "                 4. Updating Field Manager."  << G4endl;  
  fFieldManager->SetChordFinder( fChordFinder );
  fFieldManager->SetDetectorField(fMagneticField );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F01FieldSetup::~F01FieldSetup()
{
  delete fMagneticField;
  delete fChordFinder;
  delete fStepper;
  delete fFieldMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F01FieldSetup::CreateStepperAndChordFinder()
{
  delete fChordFinder;
  fChordFinder= nullptr;
   
  // Update field
  G4cout << " F01FieldSetup::CreateStepperAndChordFinder() called. " << G4endl
         << "                 1. Creating Stepper."  << G4endl;

  SetStepper();
  G4cout<<"The minimal step is equal to "<<fMinStep/mm<<" mm"<<G4endl;

  G4cout  << "                 2. Creating ChordFinder."  << G4endl;
  fChordFinder = new G4ChordFinder( fMagneticField, fMinStep,fStepper );

  G4cout  << "                 3. Updating Field Manager."  << G4endl;  
  fFieldManager->SetChordFinder( fChordFinder );
  fFieldManager->SetDetectorField(fMagneticField );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F01FieldSetup::SetStepper()
{
// Set stepper according to the stepper type

  if (fStepper) delete fStepper;

  switch ( fStepperType )
  {
     //  The new default in G4 and here ( since G4 10.4 Dec 2017 )
    case 17:      
    case 457:
    case 745:
      fStepper = new G4DormandPrince745( fEquation );
      G4cout<<"G4DormandPrince745 Stepper is chosen"<<G4endl;
      break;
     
    case 0:
      fStepper = new G4ExplicitEuler( fEquation );
      G4cout<<"G4ExplicitEuler is chosen."<<G4endl;
      break;
    case 1:
      fStepper = new G4ImplicitEuler( fEquation );
      G4cout<<"G4ImplicitEuler is chosen"<<G4endl;
      break;
    case 2:
      fStepper = new G4SimpleRunge( fEquation );
      G4cout<<"G4SimpleRunge is chosen"<<G4endl;
      break;
    case 3:
      fStepper = new G4SimpleHeum( fEquation );
      G4cout<<"G4SimpleHeum is chosen"<<G4endl;
      break;
    case 4:
      fStepper = new G4ClassicalRK4( fEquation );
      G4cout<<"G4ClassicalRK4 (default) is chosen"<<G4endl;
      break;
    case 5:
      fStepper = new G4HelixExplicitEuler( fEquation );
      G4cout<<"G4HelixExplicitEuler is chosen"<<G4endl;
      break;
    case 6:
      fStepper = new G4HelixImplicitEuler( fEquation );
      G4cout<<"G4HelixImplicitEuler is chosen"<<G4endl;
      break;
    case 7:
      fStepper = new G4HelixSimpleRunge( fEquation );
      G4cout<<"G4HelixSimpleRunge is chosen"<<G4endl;
      break;
    case 8:
      fStepper = new G4CashKarpRKF45( fEquation );
      G4cout<<"G4CashKarpRKF45 is chosen"<<G4endl;
      break;
    case 9:
      fStepper = new G4RKG3_Stepper( fEquation );
      G4cout<<"G4RKG3_Stepper is chosen"<<G4endl;
      break;
    case 10: 
       fStepper = new G4ExactHelixStepper( fEquation );   
       G4cout<<"G4ExactHelixStepper is chosen"<<G4endl;
       break;
    case 11: 
       fStepper = new G4HelixMixedStepper( fEquation );  
       G4cout<<"G4HelixMixedStepper is chosen"<<G4endl;
       break;
    case 12: 
       fStepper = new G4ConstRK4( fEquation ); 
       G4cout<<"G4ConstRK4 Stepper is chosen"<<G4endl;
       break;
    case 13:
    case 40:
      fStepper = new G4NystromRK4( fEquation );
      G4cout<<" G4NystromRK4 Stepper is chosen"<<G4endl;
      break;
    case 14:      
    case 23:
      fStepper = new G4BogackiShampine23( fEquation );
      G4cout<<"G4BogackiShampine23 Stepper is chosen"<<G4endl;
      break;

      // Other optimised 4/5th order embedded steppers
    case 15:
    case 45:       
      fStepper = new G4BogackiShampine45( fEquation );
      G4cout<<"G4BogackiShampine45 Stepper is chosen"<<G4endl;
      break;

    // case 145:
    case kTsitouras45:     
      fStepper = new G4TsitourasRK45( fEquation );
      G4cout<<"G4TsitourasRK45 Stepper is chosen"<<G4endl;
      break;      

      // Higher order embedded steppers - for very smooth fields
    case 56:
      fStepper = new G4DormandPrinceRK56( fEquation );
      G4cout<<"G4DormandPrinceRK56 Stepper is chosen"<<G4endl;
      break;
    case 78:
      fStepper = new G4DormandPrinceRK78( fEquation );
      G4cout<<"G4DormandPrinceRK78 Stepper is chosen"<<G4endl;
      break;

    default:
      // fStepper = new G4ClassicalRK4( fEquation );
      // G4cout<<"G4ClassicalRK4 Stepper (default) is chosen"<<G4endl;
      fStepper = new G4DormandPrince745( fEquation );
      G4cout<<"G4DormandPrince745 (default) Stepper is chosen"<<G4endl;
      break;      
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4VIntegrationDriver.hh"
#include "G4FSALIntegrationDriver.hh"
#include "G4RK547FEq1.hh"
#include "G4RK547FEq2.hh"
#include "G4RK547FEq3.hh"

G4VIntegrationDriver*
F01FieldSetup::CreateFSALStepperAndDriver()
{
   // using FsalStepperType = G4RK547FEq1;
  const char *methodName= "F01FieldSetup::CreateFSALStepperAndDriver()";
  if (fStepper) delete fStepper;
  fStepper = nullptr;
  
  G4cout << " F01FieldSetup::CreateFSALStepperAndDriver() called. " << G4endl;   
  G4cout << "                 1. Creating Stepper."  << G4endl;
  // auto fsalStepper = new FsalStepperType( fEquation );
  G4RK547FEq1* stepper1 = nullptr;
  G4RK547FEq2* stepper2 = nullptr;
  G4RK547FEq3* stepper3 = nullptr;

  G4cout  << "                2. Creating FSAL Driver."  << G4endl;
  G4VIntegrationDriver* fsalDriver = nullptr;
  switch ( fStepperType )
  {
    case   1:
    case 101:
       stepper1 = new G4RK547FEq1( fEquation );
       fsalDriver = new G4FSALIntegrationDriver<G4RK547FEq1>( fMinStep, stepper1 );
       G4cout  << " Stepper type '1' is G4RK547FEq1 stepper (in FSAL mode) with FSAL driver. "
               << G4endl;
       fStepper = stepper1;
       stepper1 = nullptr;
       break;
     
    case   2:
    case 102:
       stepper2= new G4RK547FEq2( fEquation );
       fsalDriver = new G4FSALIntegrationDriver<G4RK547FEq2>( fMinStep, stepper2 );
       G4cout  << " Stepper type '2' is G4RK547FEq2 stepper (in FSAL mode) with FSAL driver. "
               << G4endl;
       fStepper = stepper2;
       stepper2 = nullptr;
       break;
       
    case   3:
    case 103:
       stepper3 = new G4RK547FEq3( fEquation );       
       fsalDriver = new G4FSALIntegrationDriver<G4RK547FEq3>( fMinStep, stepper3 );
       G4cout  << " Stepper type '3' is G4RK547FEq3 stepper (in FSAL mode) with FSAL driver. "
               << G4endl;       
       fStepper = stepper3;
       stepper3 = nullptr;
       break;

    default:
       G4cout << " Warning from " << methodName << " :  stepperType (= "
              << fStepperType << " ) is unknown. " << G4endl
              << " Using value '1' instead - i.e. G4RK547FEq1 stepper. "
              << G4endl;
       stepper1 = new G4RK547FEq1( fEquation );       
       fsalDriver = new G4FSALIntegrationDriver<G4RK547FEq1>( fMinStep, stepper1 );
       fStepper = stepper1;
       stepper1 = nullptr;
       break;
  }

  delete stepper1;    stepper1 = nullptr;
  delete stepper2;    stepper2 = nullptr;              
  delete stepper3;    stepper3 = nullptr;
  
  if( fsalDriver )
     fStepper = fsalDriver->GetStepper();
  
  return fsalDriver;
}

void F01FieldSetup::CreateFSALStepperAndChordFinder()
{
  // using FsalStepperType = G4DormandPrince745; // eventually ? 
  delete fChordFinder;
  fChordFinder= nullptr;
  
  G4cout << " F01FieldSetup::CreateFSALStepperAndChordFinder() called. " << G4endl;

  auto FSALdriver= CreateFSALStepperAndDriver();
  fDriver = FSALdriver;
  G4cout<<"The minimal step is equal to "<<fMinStep/mm<<" mm"<<G4endl;

  G4cout  << "                 3. Creating ChordFinder."  << G4endl;
  fChordFinder = new G4ChordFinder( FSALdriver );  // ( fMagneticField, fMinStep, fStepper );
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F01FieldSetup::SetFieldZValue(G4double fieldStrength)
{
  // Set the value of the Global Field to fieldValue along Z

  SetFieldValue(G4ThreeVector(0, 0, fieldStrength));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F01FieldSetup::SetFieldValue(G4ThreeVector fieldVector)
{
  // Set the value of the Global Field

  if (fMagneticField) delete fMagneticField;
 
#ifdef G4VERBOSE
     G4cout << "Setting Field strength to "
            << fieldVector / gauss  << " Gauss." << G4endl;
#endif

  if (fieldVector != G4ThreeVector(0.,0.,0.))
  {
    fMagneticField = new G4UniformMagField(fieldVector);
  }
  else
  {
#ifdef G4VERBOSE
     G4cout << " Magnetic field pointer is null." << G4endl;
#endif
    // If the new field's value is Zero, signal it as below
    // so that it is not used for propagation.
    fMagneticField = 0;
  }

  // Set this as the field of the global Field Manager
  GetGlobalFieldManager()->SetDetectorField(fMagneticField);

  // Now notify equation of new field
  fEquation->SetFieldObj( fMagneticField );

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4FieldManager* F01FieldSetup::GetGlobalFieldManager()
{
  //  Utility method

  return G4TransportationManager::GetTransportationManager()
           ->GetFieldManager();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
