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
//---------------------------------------------------------------------------
//
// ClassName:   G4NeutronKiller
//
// Description: The process to kill neutrons to save CPU
//
// Author:      V.Ivanchenko 26/09/00 for HARP software
//
//----------------------------------------------------------------------------
//
// Class description:
//
// G4NeutronKiller allows to remove unwanted neutrons from simulation in 
// order to improve CPU performance. There are two parameters:
//                 
// 1) low energy threshold for neutron kinetic energy (default 0)
// 2) time limit for neutron track (default DBL_MAX) 
// 
// These parameters can be changed via Set methods or by UI commands:
//   /physics_engine/neutron/energyCut
//   /physics_engine/neutron/timeLimit
//
// If a neutron track is killed no energy deposition is added to the step 
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef NeutronKiller_h
#define NeutronKiller_h 1

#include "globals.hh"
#include "G4VDiscreteProcess.hh"
#include "G4ParticleDefinition.hh"
#include "G4Step.hh"
#include "G4Track.hh"

class G4NeutronKillerMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4NeutronKiller : public G4VDiscreteProcess
{
public:

  G4NeutronKiller(const G4String& processName = "nKiller",
		  G4ProcessType   aType =  fGeneral );

  virtual ~G4NeutronKiller();

  virtual G4bool IsApplicable(const G4ParticleDefinition&);

  void  SetTimeLimit(G4double);

  void  SetKinEnergyLimit(G4double);

  virtual G4double PostStepGetPhysicalInteractionLength( const G4Track& track,
							 G4double previousStepSize,
							 G4ForceCondition* condition);

  virtual G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);

  virtual G4double GetMeanFreePath(const G4Track&, G4double,G4ForceCondition*);

private:

  // hide assignment operator as private
  G4NeutronKiller(const G4NeutronKiller&);
  G4NeutronKiller& operator = (const G4NeutronKiller &right);

  G4double kinEnergyThreshold;
  G4double timeThreshold;
     
  G4NeutronKillerMessenger* pMess;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

