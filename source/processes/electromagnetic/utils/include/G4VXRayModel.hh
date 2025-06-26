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
//
// GEANT4 Class header file
//
//
// File name:     G4VXRayModel
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 28.04.2025
//
//
// Class Description:
//
// Abstract interface to a X-Ray production model

// -------------------------------------------------------------------
//

#ifndef G4VXRayModel_h
#define G4VXRayModel_h 1

#include "globals.hh"
#include <vector>

class G4LogicalVolume;
class G4ParticleDefinition;
class G4DynamicParticle;
class G4LossTableManager;
class G4Track;
class G4Step;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4VXRayModel
{

public:

  explicit G4VXRayModel(const G4String& nam);

  G4VXRayModel(const G4VXRayModel&);

  virtual ~G4VXRayModel();

  // return minimal beta for Cerenkov in these volumes
  G4double Initialise(std::vector<const G4LogicalVolume*>*);
  virtual G4double InitialiseModel() = 0;

  // check applicability and propose step limit, which may be DBL_MAX
  // if these methods return "false", then sampling of X-rays not possible 
  G4bool StepLimit(const G4LogicalVolume*, const G4Track&,
		   G4double preStepBeta, G4double& limit);
  virtual G4bool StepLimitForVolume(G4double& limit) = 0;
  
  // sampling is called if StepLimit(..) returns "true"
  // produced X-rays are inside vector out
  // each photon has time, position, and other parameters  within the step
  virtual void SampleXRays(std::vector<G4Track*>& out, const G4Step&) = 0;

  // for automatic documentation
  virtual void ModelDescription(std::ostream& outFile) const;

  const G4String& GetName() const { return pName; };

  G4int GetType() const { return pTypeInt; };

  //  hide assignment operator
  G4VXRayModel& operator=(const G4VXRayModel& right) = delete;

private:

  void Register();
  
protected:

  std::vector<const G4LogicalVolume*>* pLogicalVolumes{nullptr};
  G4LossTableManager* pEmManager{nullptr};
  const G4LogicalVolume* pCurrentLV{nullptr};
  const G4Track* pCurrentTrack{nullptr};
  G4double pBetaMin{1.0};
  G4double pPreStepBeta{0.0};
  G4double pMaxBetaChange{0.1};
  
  G4int pMaxPhotons{100};
  G4int pNumPhotons{0};

  G4int pTypeInt{0};
  G4int pVerbose{1};
  G4bool isMaster{false};
  const G4String pName;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

#endif
