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
// Author:        V.Ivanchenko (Vladimir.Ivantchenko@cern.ch)
// 
// Creation date: 27 September 2001
//
// Modifications: 
//
// Class Description: 
//
// Bremsstrahlung process based on the model developed  
// by Alessandra Forti, 1999, and Veronique Lefebure, 2000 
//
// -------------------------------------------------------------------

// Class description:
// Low Energy electromagnetic process,
// Further documentation available from http://www.ge.infn.it/geant4/lowE
//
// Class Description: End 

// --------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef G4LowEnergyBremsstrahlung_h
#define G4LowEnergyBremsstrahlung_h 1

#include "G4eLowEnergyLoss.hh"
#include "G4Electron.hh"
#include "G4VEMDataSet.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4Track;
class G4Step;
class G4ParticleDefinition;
class G4VParticleChange;
class G4VDataSetAlgorithm;
class G4ParticleChange;
class G4VEnergySpectrum;
class G4VCrossSectionHandler;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4LowEnergyBremsstrahlung : public G4eLowEnergyLoss

{ 
public:
 
  G4LowEnergyBremsstrahlung(const G4String& processName = "eLowEnergyBrem");
  
  ~G4LowEnergyBremsstrahlung();
  
  G4bool IsApplicable(const G4ParticleDefinition&);
  
  void PrintInfoDefinition();
  
  void BuildPhysicsTable(const G4ParticleDefinition& ParticleType);
  
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
  
  void BuildLossTable(const G4ParticleDefinition& ParticleType);
  
private:

  G4VCrossSectionHandler* crossSectionHandler;
  G4VEMDataSet* theMeanFreePath;
  G4VEnergySpectrum* theBR;

  // lower limit for generation of gamma in this model
  G4DataVector cutForSecondaryPhotons;

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
   G4int index = (track.GetMaterial())->GetIndex();
   const G4VEMDataSet* data = theMeanFreePath->GetComponent(index);
   G4double meanFreePath = data->FindValue(track.GetKineticEnergy());
   return meanFreePath; 
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  
#endif
 










