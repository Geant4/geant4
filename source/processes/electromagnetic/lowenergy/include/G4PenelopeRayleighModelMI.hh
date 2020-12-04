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
// $Id: G4PenelopeRayleighModelMI.hh 75573 2013-11-04 11:48:15Z gcosmo $
//
// Author: Luciano Pandola and Gianfranco Paternò
//
// -------------------------------------------------------------------
// History:
// 03 Dec 2009   L. Pandola   1st implementation 
// 25 May 2011   L. Pandola   Renamed (make v2008 as default Penelope)
// 27 Sep 2013   L. Pandola   Migration to MT paradigm
// 20 Aug 2017 	 G. Paternò   Molecular Interference implementation
// 24 Mar 2019 	 G. Paternò   Improved Molecular Interference implementation
// 20 Jun 2020   G. Paternò   Read qext separately and leave original atomic form factors
// 27 Aug 2020   G. Paternò   Further improvement of MI implementation
//
// -------------------------------------------------------------------
// Class description:
// Low Energy Electromagnetic Physics, Rayleigh Scattering
// with the model from Penelope, version 2008
// extended for Molecular Interference Effects
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef G4PenelopeRayleighModelMI_HH
#define G4PenelopeRayleighModelMI_HH 1

#include "globals.hh"
#include "G4VEmModel.hh"
#include "G4DataVector.hh"
#include "G4ParticleChangeForGamma.hh"

#include "G4ExtendedMaterial.hh"
#include "G4MIData.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4ParticleDefinition;
class G4DynamicParticle;
class G4MaterialCutsCouple;
class G4Material;
class G4PhysicsFreeVector;
class G4PenelopeSamplingData;

class G4PenelopeRayleighModelMI : public G4VEmModel 
{

public:  
  G4PenelopeRayleighModelMI(const G4ParticleDefinition* p = nullptr,
			    const G4String& processName = "PenRayleighMI");

  virtual ~G4PenelopeRayleighModelMI();
  
  virtual void Initialise(const G4ParticleDefinition*, const G4DataVector&) override;
  virtual void InitialiseLocal(const G4ParticleDefinition*,
			       G4VEmModel *masterModel) override;
  
  virtual G4double ComputeCrossSectionPerAtom(const G4ParticleDefinition*,
					      G4double kinEnergy,
					      G4double Z,
					      G4double A = 0,
					      G4double cut = 0,
					      G4double emax = DBL_MAX) override;
  
  //Overriding of parent's (G4VEmModel) method                                           
  virtual G4double CrossSectionPerVolume(const G4Material*,
					 const G4ParticleDefinition*,
					 G4double kineticEnergy,
					 G4double cutEnergy = 0.,
					 G4double maxEnergy = DBL_MAX) override;	
  
  virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
				 const G4MaterialCutsCouple*,
				 const G4DynamicParticle*,
				 G4double tmin,
				 G4double maxEnergy) override;
  
  void SetVerbosityLevel(G4int lev) {verboseLevel = lev;};
  G4int GetVerbosityLevel() {return verboseLevel;};
  
  //Testing purposes
  void DumpFormFactorTable(const G4Material*);
  
  //Settings
  void SetMIActive(G4bool val){fIsMIActive = val;};
  G4bool IsMIActive(){return fIsMIActive;};

private:
  G4PenelopeRayleighModelMI& operator=(const G4PenelopeRayleighModelMI &right);
  G4PenelopeRayleighModelMI(const G4PenelopeRayleighModelMI&);
  void SetParticle(const G4ParticleDefinition*);

  //Helper methods
  void ReadDataFile(G4int);   
  void ClearTables();
  void BuildFormFactorTable(const G4Material*);
  void GetPMaxTable(const G4Material*); 
  G4double GetFSquared(const G4Material*,const G4double);
  void InitializeSamplingAlgorithm(const G4Material*);
  void ReadMolInterferenceData(const G4String&,const G4String& filename="NULL");
  G4MIData* GetMIData(const G4Material*);
  void CalculateThetaAndAngFun();  
  G4double CalculateQSquared(G4double angle, G4double energy);
  G4double IntegrateFun(G4double y[], G4int n, G4double dTheta);
  void LoadKnownMIFFMaterials();


  /// Data members
  G4ParticleChangeForGamma* fParticleChange;
  const G4ParticleDefinition* fParticle;
  
  //Intrinsic energy limits of the model: cannot be extended by the parent process
  G4double fIntrinsicLowEnergyLimit;
  G4double fIntrinsicHighEnergyLimit;
  
  G4int verboseLevel;
  G4bool isInitialised;
  
  //Internal tables and manager methods
  std::map<G4int,G4PhysicsFreeVector*> *logAtomicCrossSection;
  std::map<G4int,G4PhysicsFreeVector*> *atomicFormFactor;
  std::map<G4String,G4PhysicsFreeVector*> *MolInterferenceData; //G. Paternò  
  
  G4DataVector logQSquareGrid; //log(Q^2) grid for interpolation
  std::map<const G4Material*,G4PhysicsFreeVector*> *logFormFactorTable; //log(Q^2) vs. log(F^2)
  
  G4DataVector logEnergyGridPMax; //energy grid for PMax (and originally for the x-section)
  std::map<const G4Material*,G4PhysicsFreeVector*> *pMaxTable; //E vs. Pmax
  std::map<const G4Material*,G4PenelopeSamplingData*> *samplingTable;
    
  //Used only for G4EmCalculator and Unit Tests
  G4bool fLocalTable;  
  static const G4int Ntheta = 31415;
  G4double fDTheta = {0.0001};
  G4bool fIsMIActive;	        
  G4PhysicsFreeVector* angularFunction;
  std::map<G4String,G4String> *fKnownMaterials;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif

