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
// $Id: G4EmManager.hh 66885 2013-01-16 17:37:13Z gunter $
//
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4EmManager
//
// Author:        Vladimir Ivanchenko for migration to MT
//                  
//
// Creation date: 17.05.2013
//
// Modifications:
//
//
// Class Description:
//
// A utility static class, responsable for keeping pointers
// of all shared objects and main static parameters
// for EM physics models.
//
// It is a singleton which is initialized by the master thread.
// It keeps pointers to shared objects and is read-only in run time 
//
// EM processes have to register their tables with this
// class. The responsibility of creating and deleting the tables
// remains with the energy loss classes.
//
// -------------------------------------------------------------------
//

#ifndef G4EmManager_h
#define G4EmManager_h 1

#include <map>
#include <vector>
#include "globals.hh"

class G4PhysicsTable;
class G4MaterialCutsCouple;
class G4EmManagerMessenger;
class G4ParticleDefinition;
class G4VMultipleScattering;
class G4VEmProcess;
class G4VEnergyLossProcess;
class G4EmSaturation;
class G4EmConfigurator;
class G4ElectronIonPair;
class G4LossTableBuilder;
class G4Region;

class G4EmManager
{

public:

  static G4EmManager* Instance();

  ~G4EmManager();

  //-------------------------------------------------
  // Initialisation 
  //-------------------------------------------------
  
  void PreparePhysicsTable(const G4ParticleDefinition* aParticle,
			   G4VEnergyLossProcess* p);

  void PreparePhysicsTable(const G4ParticleDefinition* aParticle,
			   G4VEmProcess* p);

  void PreparePhysicsTable(const G4ParticleDefinition* aParticle,
			   G4VMultipleScattering* p);

  void BuildPhysicsTable(const G4ParticleDefinition* aParticle);

  void BuildPhysicsTable(const G4ParticleDefinition* aParticle,
			 G4VEnergyLossProcess* p);

  void InitialiseProcess(const G4ParticleDefinition* aParticle,
			 G4VEnergyLossProcess*);

  void InitialiseProcess(const G4ParticleDefinition* aParticle,
			 G4VEmProcess*);

  void InitialiseProcess(const G4ParticleDefinition* aParticle,
			 G4VMultipleScattering*);

  //-------------------------------------------------
  // Registration of processes and tables
  //-------------------------------------------------
  
  void Register(G4VEnergyLossProcess* p);

  void Register(G4VMultipleScattering* p);

  void Register(G4VEmProcess* p);

  void RegisterExtraParticle(const G4ParticleDefinition* aParticle, 
			     G4VEnergyLossProcess* p);

  void Register(G4PhysicsTable* p);

  void DeRegister(G4PhysicsTable* p);

  //-------------------------------------------------
  // Parameter set methods
  //-------------------------------------------------
  
  void SetLossFluctuations(G4bool val);

  void SetSubCutoff(G4bool val, const G4Region* r=0);

  void SetIntegral(G4bool val);

  void SetRandomStep(G4bool val);

  void SetMinSubRange(G4double val);

  void SetMinEnergy(G4double val);

  void SetMaxEnergy(G4double val);

  void SetMaxEnergyForCSDARange(G4double val);

  void SetMaxEnergyForMuons(G4double val);

  void SetDEDXBinning(G4int val);

  void SetDEDXBinningForCSDARange(G4int val);

  void SetLambdaBinning(G4int val);

  G4int GetNumberOfBinsPerDecade() const;

  void SetStepFunction(G4double v1, G4double v2);

  void SetBuildCSDARange(G4bool val);

  void SetLPMFlag(G4bool val);

  void SetSplineFlag(G4bool val);

  void SetLinearLossLimit(G4double val);

  void SetBremsstrahlungTh(G4double val);

  void SetFactorForAngleLimit(G4double val);

  void SetVerbose(G4int val);

  //-------------------------------------------------
  // Access methods
  //-------------------------------------------------

  G4bool BuildCSDARange() const;

  G4bool LPMFlag() const;

  G4bool SplineFlag() const;

  G4double BremsstrahlungTh() const;

  G4double FactorForAngleLimit() const;

  G4double MinKinEnergy() const;

  G4double MaxKinEnergy() const;

  G4EmSaturation* EmSaturation();

  G4EmConfigurator* EmConfigurator();

  G4ElectronIonPair* ElectronIonPair();

  G4LossTableBuilder* GetTableBuilder();

private:

  //-------------------------------------------------
  // Private methods and members
  //-------------------------------------------------

  G4EmManager();

  G4VEnergyLossProcess* BuildTables(const G4ParticleDefinition* aParticle);

  void CopyTables(const G4ParticleDefinition* aParticle, 
		  G4VEnergyLossProcess*);

  void CopyDEDXTables();

  void SetParameters(const G4ParticleDefinition* aParticle,
		     G4VEnergyLossProcess*);

  G4EmManager(G4EmManager &);
  G4EmManager & operator=(const G4EmManager &right);

  static G4EmManager* theInstance;

  typedef const G4ParticleDefinition* PD;

  std::map<PD,G4VEnergyLossProcess*,std::less<PD> > loss_map;

  std::vector<G4VEnergyLossProcess*> loss_vector;
  std::vector<PD> part_vector;
  std::vector<PD> base_part_vector;
  std::vector<G4bool> tables_are_built;
  std::vector<G4bool> isActive;
  std::vector<G4PhysicsTable*> dedx_vector;
  std::vector<G4PhysicsTable*> range_vector;
  std::vector<G4PhysicsTable*> inv_range_vector;
  std::vector<G4PhysicsTable*> xsection_vector;
  std::vector<G4VMultipleScattering*> msc_vector;
  std::vector<G4VEmProcess*> emp_vector;

  G4int n_loss;
  G4int run;

  PD     firstParticle;
  PD     theElectron;

  G4bool all_tables_are_built;
  G4bool startInitialisation;

  G4bool lossFluctuationFlag;
  G4bool subCutoffFlag;
  G4bool rndmStepFlag;
  G4bool integral;
  G4bool integralActive;
  G4bool buildCSDARange;
  G4bool minEnergyActive;
  G4bool maxEnergyActive;
  G4bool maxEnergyForMuonsActive;
  G4bool stepFunctionActive;
  G4bool flagLPM;
  G4bool splineFlag;

  G4double minSubRange;
  G4double maxRangeVariation;
  G4double maxFinalStep;
  G4double minKinEnergy;
  G4double maxKinEnergy;
  G4double maxKinEnergyForMuons;
  G4double bremsTh;
  G4double factorForAngleLimit;

  G4LossTableBuilder*   tableBuilder;
  G4EmManagerMessenger* theMessenger;
  G4EmSaturation*       emSaturation;
  G4EmConfigurator*     emConfigurator;
  G4ElectronIonPair*    emElectronIonPair;

  G4int nbinsLambda;
  G4int nbinsPerDecade;
  G4int verbose;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

