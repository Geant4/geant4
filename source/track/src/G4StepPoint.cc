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
// $Id: G4StepPoint.cc 68795 2013-04-05 13:24:46Z gcosmo $
//
//
//---------------------------------------------------------------
//
//  G4StepPoint.cc
//
//  Description:
//    This class represents information associated with the
//    each end of a Step like the space/time data of the
//    particle.
//
// Contact:
//   Questions and comments to this code should be sent to
//     Hisaya Kurashige  
//
// ---------------------------------------------------------------

#include "G4StepPoint.hh"

//////////////////////////
G4StepPoint::G4StepPoint()
//////////////////////////
  : fGlobalTime(0.), fLocalTime(0.), fProperTime(0.),
    fKineticEnergy(0.), fVelocity(0.), fpTouchable(0), fpMaterial(0),
    fpMaterialCutsCouple(0), fpSensitiveDetector(0), fSafety(0.),
    fStepStatus(fUndefined),fpProcessDefinedStep(0), 
    fMass(0.), fCharge(0.), fMagneticMoment(0.),
    fWeight(0.)
{
}

//////////////////////////
G4StepPoint::G4StepPoint(const G4StepPoint &right) :
//////////////////////////
  fPosition(right.fPosition),
  fGlobalTime(right.fGlobalTime),
  fLocalTime(right.fLocalTime),
  fProperTime(right.fProperTime),
  fMomentumDirection(right.fMomentumDirection),
  fKineticEnergy(right.fKineticEnergy),
  fVelocity(right.fVelocity),
  fpTouchable(right.fpTouchable),
  fpMaterial(right.fpMaterial),
  fpMaterialCutsCouple(right.fpMaterialCutsCouple),
  fpSensitiveDetector(right.fpSensitiveDetector),
  fSafety(right.fSafety),
  fPolarization(right.fPolarization),
  fStepStatus(right.fStepStatus),
  fpProcessDefinedStep(right.fpProcessDefinedStep),
  fMass(right.fMass),
  fCharge(right.fCharge),
  fMagneticMoment(right.fMagneticMoment),
  fWeight(right.fWeight)
{
}

//////////////////////////
G4StepPoint & G4StepPoint::operator=(const G4StepPoint &right)
{
  if (this != &right) {
    fPosition     = right.fPosition;
    fGlobalTime   = right.fGlobalTime;
    fLocalTime    = right.fLocalTime;
    fProperTime   = right.fProperTime;
    fMomentumDirection = right.fMomentumDirection;
    fKineticEnergy = right.fKineticEnergy;
    fVelocity     = right.fVelocity;
    fpTouchable   = right.fpTouchable;
    fpMaterial    = right.fpMaterial;
    fpMaterialCutsCouple = right.fpMaterialCutsCouple;
    fpSensitiveDetector  = right.fpSensitiveDetector;
    fSafety       = right.fSafety;
    fPolarization = right.fPolarization;
    fStepStatus   = right.fStepStatus;
    fpProcessDefinedStep = right.fpProcessDefinedStep;
    fMass         = right.fMass;
    fCharge       = right.fCharge;
    fMagneticMoment = right.fMagneticMoment;
    fWeight       = right.fWeight;
  }
  return *this;
}
