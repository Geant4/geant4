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
// $Id: G4NuclearStopping.hh,v 1.2 2009-11-10 19:25:47 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -----------------------------------------------------------------------------
//
// GEANT4 Class header file
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 20 July 2009
// 
// Modified:
//
//------------------------------------------------------------------------------
//

// class description
//
//  The class simulates nuclear stopping due to multiple scattering
//
// class description - end

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef G4NuclearStopping_h
#define G4NuclearStopping_h 1

#include "G4VEmProcess.hh"
#include "globals.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ParticleDefinition.hh"
#include "G4GPILSelection.hh"
#include "G4ParticleChangeForLoss.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4ICRU49NuclearStoppingModel;

class G4NuclearStopping : public G4VEmProcess
{

public:    // with description

  G4NuclearStopping(const G4String& processName="nuclearStopping");

  virtual ~G4NuclearStopping();

  // returns true for charged particles, false otherwise
  G4bool IsApplicable (const G4ParticleDefinition& p);

  // implementation of pure virtual method
  G4double AlongStepGetPhysicalInteractionLength(const G4Track& track,
						 G4double  previousStepSize,
						 G4double  currentMinimumStep,
						 G4double& proposedSafety,
						 G4GPILSelection* selection);

  // implementation of energy loss along step
  G4VParticleChange* AlongStepDoIt(const G4Track& track,
				   const G4Step&  step);

  // Print few lines of informations about the process: validity range,
  void PrintInfo();

protected:

  // This function initialise process
  void InitialiseProcess(const G4ParticleDefinition*);

private:       

  G4ParticleChangeForLoss nParticleChange;

  G4ICRU49NuclearStoppingModel* modelICRU49;

  G4bool   isInitialized;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
