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
// $Id: G4StepPoint.cc,v 1.8 2002-12-04 22:46:00 tsasaki Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
//     Katsuya Amako  (e-mail: Katsuya.Amako@kek.jp)
//     Takashi Sasaki (e-mail: Takashi.Sasaki@kek.jp)
//
// ---------------------------------------------------------------

#include "G4StepPoint.hh"

//////////////////////////
G4StepPoint::G4StepPoint():
//////////////////////////
  fpTouchable(0),fpMaterial(0),fpProcessDefinedStep(0),fCharge(0.)
{

}

//////////////////////////
G4StepPoint::G4StepPoint(const G4StepPoint &right):
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
  fSafety(right.fSafety),
  fPolarization(right.fPolarization),
  fStepStatus(right.fStepStatus),
  fpProcessDefinedStep(right.fpProcessDefinedStep),
  fMass(right.fMass),
  fCharge(right.fCharge),
  fWeight(right.fWeight)
{}


//////////////////////////
G4StepPoint & G4StepPoint::operator=(const G4StepPoint &right)
//////////////////////////
{
  if (this != &right) {
    fPosition = right.fPosition;
    fGlobalTime = right.fGlobalTime;
    fLocalTime = right.fLocalTime;
    fProperTime = right.fProperTime;
    fMomentumDirection = right.fMomentumDirection;
    fKineticEnergy = right.fKineticEnergy;
    fpTouchable = right.fpTouchable;
    fpMaterial = right.fpMaterial;
    fSafety = right.fSafety;
    fPolarization = right.fPolarization;
    fStepStatus = right.fStepStatus;
    fpProcessDefinedStep = right.fpProcessDefinedStep;
    fMass = right.fMass;
    fCharge = right.fCharge;
    fWeight = right.fWeight;
    fVelocity = right.fVelocity;
  }
  return *this;
}
