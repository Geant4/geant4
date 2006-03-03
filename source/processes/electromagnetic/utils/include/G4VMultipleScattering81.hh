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
// $Id: G4VMultipleScattering81.hh,v 1.2 2006-03-03 14:11:45 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4VMultipleScattering81
//
// Author:        Vladimir Ivanchenko 
//
// Creation date: 01.03.2006
//
// Modifications:
//
//
//
// Class Description:
//
// It is the generic process of multiple scattering it includes common
// part of calculations for all charged particles
//
//

// -------------------------------------------------------------------
//

#ifndef G4VMultipleScattering81_h
#define G4VMultipleScattering81_h 1

#include "G4VMultipleScattering.hh"
#include "globals.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4VEmModel.hh"
#include "G4VParticleChange.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4VMultipleScattering81 : public G4VMultipleScattering
{
public:

  G4VMultipleScattering81(const G4String& name = "msc81",
			  G4ProcessType type = fElectromagnetic):
    G4VMultipleScattering(name) {};

  virtual ~G4VMultipleScattering81() {};

  //------------------------------------------------------------------------
  // Virtual methods to be implemented for the concrete model
  //------------------------------------------------------------------------

public:

  //------------------------------------------------------------------------
  // Generic methods common to all models
  //------------------------------------------------------------------------

  // This method is used for tracking, it returns step limit
  virtual G4double GetContinuousStepLimit(const G4Track& track,
                                        G4double previousStepSize,
                                        G4double currentMinimalStep,
                                        G4double& currentSafety);

  virtual G4VParticleChange* AlongStepDoIt(const G4Track&, const G4Step&);

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VMultipleScattering81::GetContinuousStepLimit(
                                          const G4Track& track,
                                                G4double,
                                                G4double currentMinimalStep,
                                                G4double& currentSafety)
{
  DefineMaterial(track.GetMaterialCutsCouple());
  SelectModel(track.GetKineticEnergy());
  truePathLength = currentModel->ComputeTruePathLengthLimit(
		   track, theLambdaTable, currentMinimalStep, currentSafety);
  if (truePathLength < currentMinimalStep) valueGPILSelectionMSC = CandidateForSelection;
  geomPathLength = currentModel->ComputeGeomPathLength(truePathLength);
  return geomPathLength;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4VParticleChange* G4VMultipleScattering81::AlongStepDoIt(
                                                        const G4Track&,
                                                        const G4Step& step)
{
  G4double geomStepLength = step.GetStepLength();
  if(geomStepLength == geomPathLength) trueStepLength = truePathLength;
  else   trueStepLength = currentModel->TrueStepLength(geomStepLength);
  fParticleChange.ProposeTrueStepLength(trueStepLength);
  return &fParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
