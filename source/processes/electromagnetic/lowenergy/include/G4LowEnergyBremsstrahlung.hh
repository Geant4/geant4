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
// ------------------------------------------------------------
//      GEANT 4 class header file
//      CERN Geneva Switzerland
//
//      ------------ G4LowEnergyBremsstrahlung physics process ------
//                     by A.Forti  1999/03/27 19:18:13
//
// Class description:
// Low Energy electromagnetic process, Bremsstrahlung
// Further documentation available from http://www.ge.infn.it/geant4/lowE
//
// Class Description: End 
//
// 18.04.2000 V.Lefebure First implementation of continuous energy loss.
// 21.08.2001 V.Ivanchenko new design iteration
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef G4LowEnergyBremsstrahlung_h
#define G4LowEnergyBremsstrahlung_h 1

#include "G4eLowEnergyLoss.hh"
#include "G4Electron.hh"
#include "G4VEMSecondaryGenerator.hh"
#include "G4VEMDataSet.hh"
#include "G4CrossSectionHandler.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4Track;
class G4Step;
class G4ParticleDefinition;
class G4VParticleChange;
class G4VEMDataSet;
class G4VDataSetAlgorithm;
class G4ParticleChange;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4LowEnergyBremsstrahlung : public G4eLowEnergyLoss

{ 
public:
 
  G4LowEnergyBremsstrahlung(const G4String& processName = "eLowEnergyBrem");
  
  ~G4LowEnergyBremsstrahlung();
  
  G4bool IsApplicable(const G4ParticleDefinition&);

  //  void SetCutForLowEnSecPhotons(G4double cut) 
  //     {if(theGenerator) theGenerator->SetMinSecondaryEnergy(cut);};
  
  void PrintInfoDefinition();
  
  void BuildPhysicsTable(const G4ParticleDefinition& ParticleType);
  
  void BuildLossTable(const G4ParticleDefinition& ParticleType);
  
  G4double GetMeanFreePath(const G4Track& track,
			         G4double previousStepSize,
			         G4ForceCondition* condition );
 
  G4VParticleChange* PostStepDoIt(const G4Track& track,         
				  const G4Step&  step);                 
  
    
private:

  // Hide copy constructor and assignment operator as private 
  G4LowEnergyBremsstrahlung(const G4LowEnergyBremsstrahlung& );
  G4LowEnergyBremsstrahlung& operator = 
                             (const G4LowEnergyBremsstrahlung& right);
    
private:

  G4CrossSectionHandler* shellCrossSectionHandler;
  G4VEMDataSet* theMeanFreePathData;
  G4VEMSecondaryGenerator* theGenerator;

  // lower limit for generation of gamma in this model
  G4DataVector cutForLowEnergySecondaryPhotons;
  G4DataVector tmax;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4bool G4LowEnergyBremsstrahlung::IsApplicable(
                            const G4ParticleDefinition& particle)
{
   return(  (&particle == G4Electron::Electron())
          /////////////||(&particle == G4Positron::Positron())
	   );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4LowEnergyBremsstrahlung::GetMeanFreePath(
                            const G4Track& track,
                                  G4double,
                                  G4ForceCondition* cond)
{
   *cond = NotForced;
   G4double meanFreePath = theMeanFreePathData->FindValue(
                           track.GetKineticEnergy(),
                           (G4int)(track.GetMaterial())->GetIndex());
   return meanFreePath; 
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4VParticleChange* G4LowEnergyBremsstrahlung::PostStepDoIt(
                                                 const G4Track& track,
                                                 const G4Step&  step)
{
  aParticleChange.Initialize(track);
  const G4Material* mat = track.GetMaterial();
  G4int Z = shellCrossSectionHandler->SelectRandomAtom(mat,
                                      track.GetKineticEnergy());

  theGenerator->GenerateSecondary(track.GetDynamicParticle(), 
     &aParticleChange, Z, -1,
     cutForLowEnergySecondaryPhotons[mat->GetIndex()], 
     DBL_MAX);

  return G4VContinuousDiscreteProcess::PostStepDoIt(track, step);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  
#endif
 










