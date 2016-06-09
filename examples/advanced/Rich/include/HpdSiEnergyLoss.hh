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
// Rich advanced example for Geant4
// HpdSiEnergyLoss.hh for Rich of LHCb
// History:
// Created: Sajan Easo (Sajan.Easo@cern.ch)
// Revision and changes: Patricia Mendez (Patricia.Mendez@cern.ch)
/////////////////////////////////////////////////////////////////////////////
#ifndef HpdSiEnergyLoss_h
#define HpdSiEnergyLoss_h 1
#include "globals.hh"
#include "G4VEnergyLoss.hh"
#include "G4VParticleChange.hh"
#include "G4DynamicParticle.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ParticleDefinition.hh"

class HpdSiEnergyLoss  : public G4VEnergyLoss 
{
public:  

  HpdSiEnergyLoss(const G4String& materialName , 
		  const G4String& processName);


  virtual ~HpdSiEnergyLoss();

  G4bool IsApplicable(const G4ParticleDefinition&);
  // true for all charged particles
 G4double GetContinuousStepLimit(const G4Track& track,
                                G4double previousStepSize,
                                G4double currentMinimumStep,
                                G4double& currentSafety);
 G4VParticleChange* AlongStepDoIt(const G4Track& aTrack,
                                     const G4Step& aStep) ;
 G4double GetMeanFreePath(const G4Track& track,
                          G4double previousStepSize,
			  G4ForceCondition* condition);
 G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
                                            const G4Step& aStep);
private:

   // hide  assignment and copy operators
     HpdSiEnergyLoss ( HpdSiEnergyLoss &); 
     HpdSiEnergyLoss & operator=(const  HpdSiEnergyLoss &right);
 
  //now the data members
    G4int fMatIndex;
    G4double MinKineticEnergy ;
    G4double MipEnergy;
    G4double finalRangeforStep;   
    G4String ElossMaterialName;

};
#endif



