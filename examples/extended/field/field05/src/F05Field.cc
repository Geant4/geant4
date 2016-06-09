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
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "F05Field.hh"

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"

#include "G4EqEMFieldWithSpin.hh"
#include "G4ChordFinder.hh"
#include "G4PropagatorInField.hh"

#include "G4MagIntegratorStepper.hh"
#include "G4ClassicalRK4.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F05Field::F05Field() : G4ElectroMagneticField()
{
  fEquation = new G4EqEMFieldWithSpin(this);

  G4FieldManager* fieldMgr
   = G4TransportationManager::GetTransportationManager()->GetFieldManager();

  fieldMgr->SetDetectorField(this);

  fStepper = new G4ClassicalRK4(fEquation,12);

  G4double minStep           = 0.01*mm;

  fChordFinder = new G4ChordFinder((G4MagneticField*)this,minStep,fStepper);
 
  // Set accuracy parameters
  G4double deltaChord        = 3.0*mm;
  fChordFinder->SetDeltaChord( deltaChord );

  G4double deltaOneStep      = 0.01*mm; 
  fieldMgr->SetAccuraciesWithDeltaOneStep(deltaOneStep);

  G4double deltaIntersection = 0.1*mm; 
  fieldMgr->SetDeltaIntersection(deltaIntersection);

  G4TransportationManager* fTransportManager =
         G4TransportationManager::GetTransportationManager();

  fieldPropagator = fTransportManager->GetPropagatorInField();

  G4double epsMin            = 2.5e-7*mm;
  G4double epsMax            = 0.05*mm;
 
  fieldPropagator->SetMinimumEpsilonStep(epsMin);
  fieldPropagator->SetMaximumEpsilonStep(epsMax);
 
  fieldMgr->SetChordFinder(fChordFinder);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F05Field::~F05Field()
{
  if (fEquation)    delete fEquation;
  if (fStepper)     delete fStepper;
  if (fChordFinder) delete fChordFinder;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F05Field::GetFieldValue( const G4double Point[3],G4double* Bfield ) const
{
  // Point[0],Point[1],Point[2] are x-, y-, z-cordinates 

  const G4double Bz = 0.24*tesla;
  const G4double Er = 2.113987E+6*volt/m;

  G4double Ex,Ey;
 
  G4double posR = std::sqrt(std::pow(Point[0],2) + std::pow(Point[1],2));
  G4double cos_theta, sin_theta;

  if (posR>0){
     cos_theta = Point[0]/(G4double)posR;
     sin_theta = Point[1]/(G4double)posR;
     Ex = -1*Er*cos_theta;//apply radial electric field
     Ey = -1*Er*sin_theta;
  }else{
     Ex=0;
     Ey=0;
  }
  
  Bfield[0]=0;
  Bfield[1]=0;
  Bfield[2]=Bz;

  Bfield[3]=Ex;
  Bfield[4]=Ey;
  Bfield[5]=0;

  return;
}
