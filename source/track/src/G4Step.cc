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
// G4Step class implementation
//
// Authors:
//   Katsuya Amako  (e-mail: Katsuya.Amako@kek.jp)
//   Takashi Sasaki (e-mail: Takashi.Sasaki@kek.jp)
// Revisions:
//   Hisaya Kurashige, 1998-2007
// --------------------------------------------------------------------

#include "G4Step.hh"

// --------------------------------------------------------------------
G4Step::G4Step()
{
  fpPreStepPoint  = new G4StepPoint();
  fpPostStepPoint = new G4StepPoint();

  secondaryInCurrentStep = new std::vector<const G4Track*>;
}

// --------------------------------------------------------------------
G4Step::~G4Step()
{
  delete fpPreStepPoint;
  delete fpPostStepPoint;

  secondaryInCurrentStep->clear();
  delete secondaryInCurrentStep;

  if(fSecondary != nullptr)
  {
    fSecondary->clear();
    delete fSecondary;
  }
}

// --------------------------------------------------------------------
G4Step::G4Step(const G4Step& right)
  : fTotalEnergyDeposit(right.fTotalEnergyDeposit)
  , fNonIonizingEnergyDeposit(right.fNonIonizingEnergyDeposit)
  , fStepLength(right.fStepLength)
  , fpTrack(right.fpTrack)
  , fpSteppingControlFlag(right.fpSteppingControlFlag)
  , fFirstStepInVolume(right.fFirstStepInVolume)
  , fLastStepInVolume(right.fLastStepInVolume)
  , nSecondaryByLastStep(right.nSecondaryByLastStep)
  , secondaryInCurrentStep(right.secondaryInCurrentStep)
  , fpVectorOfAuxiliaryPointsPointer(right.fpVectorOfAuxiliaryPointsPointer)
{
  if(right.fpPreStepPoint != nullptr)
  {
    fpPreStepPoint = new G4StepPoint(*(right.fpPreStepPoint));
  }
  else
  {
    fpPreStepPoint = new G4StepPoint();
  }
  if(right.fpPostStepPoint != nullptr)
  {
    fpPostStepPoint = new G4StepPoint(*(right.fpPostStepPoint));
  }
  else
  {
    fpPostStepPoint = new G4StepPoint();
  }

  if(right.fSecondary != nullptr)
  {
    fSecondary = new G4TrackVector(*(right.fSecondary));
  }
  else
  {
    fSecondary = new G4TrackVector();
  }

  // secondaryInCurrentStep is cleared
  secondaryInCurrentStep = new std::vector<const G4Track*>;
}

// --------------------------------------------------------------------
G4Step& G4Step::operator=(const G4Step& right)
{
  if(this != &right)
  {
    fTotalEnergyDeposit              = right.fTotalEnergyDeposit;
    fNonIonizingEnergyDeposit        = right.fNonIonizingEnergyDeposit;
    fStepLength                      = right.fStepLength;
    fpTrack                          = right.fpTrack;
    fpSteppingControlFlag            = right.fpSteppingControlFlag;
    fFirstStepInVolume               = right.fFirstStepInVolume;
    fLastStepInVolume                = right.fLastStepInVolume;
    nSecondaryByLastStep             = right.nSecondaryByLastStep;
    secondaryInCurrentStep           = right.secondaryInCurrentStep;
    fpVectorOfAuxiliaryPointsPointer = right.fpVectorOfAuxiliaryPointsPointer;

    delete fpPreStepPoint;

    if(right.fpPreStepPoint != nullptr)
    {
      fpPreStepPoint = new G4StepPoint(*(right.fpPreStepPoint));
    }
    else
    {
      fpPreStepPoint = new G4StepPoint();
    }

    delete fpPostStepPoint;

    if(right.fpPostStepPoint != nullptr)
    {
      fpPostStepPoint = new G4StepPoint(*(right.fpPostStepPoint));
    }
    else
    {
      fpPostStepPoint = new G4StepPoint();
    }

    if(fSecondary != nullptr)
    {
      fSecondary->clear();
      delete fSecondary;
    }
    if(right.fSecondary != nullptr)
    {
      fSecondary = new G4TrackVector(*(right.fSecondary));
    }
    else
    {
      fSecondary = new G4TrackVector();
    }

    // secondaryInCurrentStep is not copied
    if(secondaryInCurrentStep != nullptr)
    {
      secondaryInCurrentStep->clear();
      delete secondaryInCurrentStep;
    }
    secondaryInCurrentStep = new std::vector<const G4Track*>;
  }
  return *this;
}

// --------------------------------------------------------------------
G4ThreeVector G4Step::GetDeltaMomentum() const
{
  static G4ThreadLocal G4bool isFirstTime = true;
  if(isFirstTime)
  {
    isFirstTime = false;
#ifdef G4VERBOSE
    G4Exception("G4Step::GetDeltaMomentum()", "Warning", JustWarning,
                "This method is obsolete and will be removed soon");
#endif
  }

  return fpPostStepPoint->GetMomentum() - fpPreStepPoint->GetMomentum();
}

// --------------------------------------------------------------------
G4double G4Step::GetDeltaEnergy() const
{
  static G4ThreadLocal G4bool isFirstTime = true;
  if(isFirstTime)
  {
    isFirstTime = false;
#ifdef G4VERBOSE
    G4Exception("G4Step::GetDeltaEnergy()", "Warning", JustWarning,
                "This method is obsolete and will be removed soon");
#endif
  }

  return fpPostStepPoint->GetKineticEnergy() -
         fpPreStepPoint->GetKineticEnergy();
}

// --------------------------------------------------------------------
const std::vector<const G4Track*>* G4Step::GetSecondaryInCurrentStep() const
{
  secondaryInCurrentStep->clear();
  std::size_t nSecondary = fSecondary->size();
  for(std::size_t i = nSecondaryByLastStep; i < nSecondary; ++i)
  {
    secondaryInCurrentStep->push_back((*fSecondary)[i]);
  }
  return secondaryInCurrentStep;
}
