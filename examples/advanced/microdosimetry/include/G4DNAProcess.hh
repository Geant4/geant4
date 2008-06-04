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
// -------------------------------------------------------------------
// $Id: G4DNAProcess.hh,v 1.1 2008-06-04 12:58:24 sincerti Exp $
// -------------------------------------------------------------------

#ifndef G4DNAPROCESS_HH
#define G4DNAPROCESS_HH 1

#include "globals.hh"
#include "G4VDiscreteProcess.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ParticleDefinition.hh"
#include "G4VParticleChange.hh"
#include "G4DynamicParticle.hh"
#include "G4ThreeVector.hh"
#include "G4FinalStateProduct.hh"
#include "G4String.hh"
#include <vector> 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

template < class TCrossSection, class TFinalState, char const* name>
class G4DNAProcess : public G4VDiscreteProcess
{
public:

  G4DNAProcess(const G4String& processName = name): G4VDiscreteProcess(processName) { ; }
  
  ~G4DNAProcess() { ; }

  virtual G4bool IsApplicable(const G4ParticleDefinition&) { return true; } 
  
  virtual void BuildPhysicsTable(const G4ParticleDefinition&) { ; }
 
  virtual G4VParticleChange* PostStepDoIt(const G4Track& track, const G4Step& step);
 
  G4double DumpMeanFreePath(const G4Track& aTrack, 
			    G4double previousStepSize, 
			    G4ForceCondition* condition) 
  { return GetMeanFreePath(aTrack, previousStepSize, condition); }

protected:
  
  virtual G4double GetMeanFreePath(const G4Track& track, 
			   G4double previousStepSize, 
			   G4ForceCondition* condition);
private:

  G4DNAProcess& operator=(const G4DNAProcess& right);
  G4DNAProcess(const G4DNAProcess& );

  TCrossSection crossSection;
  TFinalState finalState;
};
#include "G4DNAProcess.icc"
#endif
