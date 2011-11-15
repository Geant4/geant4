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
// $Id: StepLimiterPerRegion.hh,v 1.3 2008-08-05 10:38:35 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef StepLimiterPerRegion_h
#define StepLimiterPerRegion_h 1

#include "globals.hh"
#include "G4VDiscreteProcess.hh"
#include "G4ParticleDefinition.hh"
#include "G4Step.hh"

class StepLimiterMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class StepLimiterPerRegion : public G4VDiscreteProcess
{
public:

  StepLimiterPerRegion(const G4String& processName = "UserMaxStep");
  virtual ~StepLimiterPerRegion();

  G4bool IsApplicable(const G4ParticleDefinition&);

  void SetMaxStep(G4double);

  G4double GetMaxStep() {return MaxChargedStep;};

  G4double PostStepGetPhysicalInteractionLength( const G4Track& track,
			                       G4double previousStepSize,
			                       G4ForceCondition* condition);

  G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);

  G4double GetMeanFreePath(const G4Track&, G4double, G4ForceCondition*)
  {return 0.;};    

private:

  StepLimiterPerRegion & operator=(const StepLimiterPerRegion &right);
  StepLimiterPerRegion(const StepLimiterPerRegion&);

  G4double MaxChargedStep;
  G4double ProposedStep;

  StepLimiterMessenger*  pMess;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

