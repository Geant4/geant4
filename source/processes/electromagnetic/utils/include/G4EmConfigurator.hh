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
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
// File name:     G4EmConfigurator
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 14.07.2008
//
// Modifications:
//
// Class Description:
//
// This class provides configuration EM models for 
// particles/processes/regions
//

// -------------------------------------------------------------------
//

#ifndef G4EmConfigurator_h
#define G4EmConfigurator_h 1

#include "globals.hh"
#include "G4VEmModel.hh"
#include "G4VEmFluctuationModel.hh"
#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4VEnergyLossProcess;
class G4VEmProcess;
class G4VMultipleScattering;
class G4TransportationWithMsc;

class G4EmConfigurator 
{
public: 
  
  explicit G4EmConfigurator(G4int verboseLevel = 1);
 
  ~G4EmConfigurator();

  // Set EM model for particle type and process to 
  // be active for the G4Region and energy interval
  // The model will be added to the list 
  //
  void SetExtraEmModel(const G4String& particleName,
                       const G4String& processName,
                       G4VEmModel*,
                       const G4String& regionName = "",
                       G4double emin = 0.0,
                       G4double emax = DBL_MAX,
                       G4VEmFluctuationModel* fm = nullptr); 

  // Add all previously declared models to corresponding processes
  // Can be called in ConstructPhysics
  //
  void AddModels();

  // These methods called by G4LossTableManager
  //
  void PrepareModels(const G4ParticleDefinition* aParticle,
                     G4VEnergyLossProcess* p);

  void PrepareModels(const G4ParticleDefinition* aParticle,
                     G4VEmProcess* p);

  void PrepareModels(const G4ParticleDefinition* aParticle,
                     G4VMultipleScattering* p,
                     G4TransportationWithMsc* trans = nullptr);

  void Clear();

  inline void SetVerbose(G4int value);

  // hide assignment operator
  G4EmConfigurator & operator=(const G4EmConfigurator &right) = delete;
  G4EmConfigurator(const G4EmConfigurator&) = delete;

private:

  void SetModelForRegion(G4VEmModel* model,
                         G4VEmFluctuationModel* fm,
                         const G4Region* reg,
                         const G4String& particleName,
                         const G4String& processName,
                         G4double emin,
                         G4double emax);

  G4bool UpdateModelEnergyRange(G4VEmModel* mod,
                                G4double emin, G4double emax);

  std::vector<G4VEmModel*> models;  
  std::vector<G4VEmFluctuationModel*> flucModels;  
  std::vector<G4String> particles;  
  std::vector<G4String> processes;  
  std::vector<G4String> regions;  
  std::vector<G4double> lowEnergy;
  std::vector<G4double> highEnergy;
  
  G4int index;
  G4int verbose;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void G4EmConfigurator::SetVerbose(G4int value)
{
  verbose = value;
}

#endif








