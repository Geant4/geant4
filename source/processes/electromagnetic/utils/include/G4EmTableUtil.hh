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
// Geant4 header G4EmTableUtil
//
// Author V.Ivanchenko 14.03.2022
//
// Utilities used at initialisation of EM physics
//

#ifndef G4EmTableUtil_h
#define G4EmTableUtil_h 1

#include "globals.hh"
#include "G4PhysicsTable.hh"
#include "G4VMultipleScattering.hh"
#include "G4VEmProcess.hh"
#include "G4VEnergyLossProcess.hh"
#include "G4EmModelManager.hh"
#include "G4LossTableBuilder.hh"
#include "G4ParticleDefinition.hh"
#include "G4MscStepLimitType.hh"
#include "G4DataVector.hh"

class G4EmTableUtil
{
public:

  static const G4DataVector*
  PrepareEmProcess(G4VEmProcess* proc,
                   const G4ParticleDefinition* part,
                   const G4ParticleDefinition* secPart,
		   G4EmModelManager* modelManager,
                   const G4double& maxKinEnergy,
		   G4int& secID, G4int& tripletID,
                   G4int& mainSec, const G4int& verb,
                   const G4bool& master);

  static void BuildEmProcess(G4VEmProcess* proc,
                             const G4VEmProcess* masterProc,
			     const G4ParticleDefinition* firstPart,
			     const G4ParticleDefinition* part,
			     const G4int nModels, const G4int verb,
                             const G4bool master, const G4bool isLocked,
                             const G4bool toBuild, G4bool& baseMat);

  static void BuildLambdaTable(G4VEmProcess* proc,
                               const G4ParticleDefinition* part,
                               G4EmModelManager* modelManager,
			       G4LossTableBuilder* bld,
                               G4PhysicsTable* theLambdaTable,
                               G4PhysicsTable* theLambdaTablePrim,
                               const G4double minKinEnergy,
                               const G4double minKinEnergyPrim,
                               const G4double maxKinEnergy,
                               const G4double scale,
                               const G4int verbose,
                               const G4bool startFromNull,
                               const G4bool splineFlag);

  static void BuildLambdaTable(G4VEnergyLossProcess* proc,
                               const G4ParticleDefinition* part,
                               G4EmModelManager* modelManager,
			       G4LossTableBuilder* bld,
                               G4PhysicsTable* theLambdaTable,
                               const G4DataVector* theCuts,
                               const G4double minKinEnergy,
                               const G4double maxKinEnergy,
                               const G4double scale,
                               const G4int verbose,
                               const G4bool splineFlag);

  static const G4ParticleDefinition* CheckIon(
                               G4VEnergyLossProcess* proc,
                               const G4ParticleDefinition* part,
                               const G4ParticleDefinition* particle,
                               const G4int verboseLevel, G4bool& isIon);

  static void UpdateModels(G4VEnergyLossProcess* proc,
			   G4EmModelManager* modelManager,
                           const G4double maxKinEnergy,
                           const G4int nModels,
                           G4int& secID, G4int& biasID,
                           G4int& mainSecondaries, const G4bool baseMat,
                           const G4bool isMaster, const G4bool useAGen);

  static void BuildLocalElossProcess(G4VEnergyLossProcess* proc,
				     const G4VEnergyLossProcess* masterProc,
				     const G4ParticleDefinition* part,
                                     const G4int nModels);

  static void BuildDEDXTable(G4VEnergyLossProcess* proc,
			     const G4ParticleDefinition* part,
			     G4EmModelManager* modelManager,
			     G4LossTableBuilder* bld,
			     G4PhysicsTable* table,
			     const G4double minKinEnergy,
			     const G4double maxKinEnergy,
			     const G4int nbins,
			     const G4int verbose,
			     const G4EmTableType tType,
			     const G4bool splineFlag);

  static void PrepareMscProcess(G4VMultipleScattering* proc,
                                const G4ParticleDefinition& part,
			        G4EmModelManager* modelManager,
			        G4MscStepLimitType& stepLimit,
                                G4double& facrange,
			        G4bool& latDisplacement, G4bool& master,
			        G4bool& isIon, G4bool& baseMat);

  static void BuildMscProcess(G4VMultipleScattering* proc,
                              const G4VMultipleScattering* masterProc,
		              const G4ParticleDefinition& part,
		              const G4ParticleDefinition* firstPart,
		              G4int nModels, G4bool master);

  static G4bool StoreMscTable(G4VMultipleScattering* proc,
                              const G4ParticleDefinition* part,
                              const G4String& directory,
			      const G4int nModels, const G4int verb,
                              const G4bool ascii);

  static G4bool StoreTable(G4VProcess*, const G4ParticleDefinition*, 
                           G4PhysicsTable*, const G4String& dir,
                           const G4String& tname, G4int verb,
                           G4bool ascii);

  static G4bool RetrieveTable(G4VProcess* ptr,
                              const G4ParticleDefinition* part, 
                              G4PhysicsTable* aTable, 
                              const G4String& dir, const G4String& tname,
                              const G4int verb, const G4bool ascii,
                              const G4bool spline);

};

#endif


