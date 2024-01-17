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
// Code developed by
// Silvia Pozzi (1), silvia.pozzi@iss.it
// Barbara Caccia (1), barbara.caccia@iss.it
// Carlo Mancini Terracciano (2), carlo.mancini.terracciano@roma1.infn.it
// (1) Istituto Superiore di Sanita' and INFN Roma, Italy
// (2) Univ. La Sapienza and INFN Roma, Italy

#ifndef StepMax_h
#define StepMax_h 1

#include "globals.hh"
#include "G4VDiscreteProcess.hh"
#include "G4ParticleDefinition.hh"
#include "G4Step.hh"

class StepMaxMessenger;

/////////////////////////////////////////////////////////////////////////////
class StepMax : public G4VDiscreteProcess
{
public:     
  
  StepMax(const G4String& processName ="UserStepMax");
  ~StepMax();
  
  G4bool   IsApplicable(const G4ParticleDefinition&) override;
  void     SetMaxStep(G4double);
  inline G4double GetMaxStep() {return MaxChargedStep;}
  
  G4double PostStepGetPhysicalInteractionLength( const G4Track& track,
						 G4double   previousStepSize,
						 G4ForceCondition* condition) override;
  
  G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&) override;
  
  G4double GetMeanFreePath(const G4Track&, G4double, G4ForceCondition*) override
  {return 0.;}     // it is not needed here !
  
private:

	 G4double    MaxChargedStep;
	 StepMaxMessenger* pMess;
};

#endif

