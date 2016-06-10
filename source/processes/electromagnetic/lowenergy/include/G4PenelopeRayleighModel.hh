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
// $Id: G4PenelopeRayleighModel.hh 75573 2013-11-04 11:48:15Z gcosmo $
//
// Author: Luciano Pandola
//
// History:
// -----------
// 03 Dec 2009   L. Pandola   1st implementation. 
// 25 May 2011   L. Pandola   Renamed (make v2008 as default Penelope)
// 27 Sep 2013   L. Pandola  Migration to MT paradigm
//
// -------------------------------------------------------------------
//
// Class description:
// Low Energy Electromagnetic Physics, Rayleigh Scattering
// with the model from Penelope, version 2008
// -------------------------------------------------------------------

#ifndef G4PENELOPERAYLEIGHMODEL_HH
#define G4PENELOPERAYLEIGHMODEL_HH 1

#include "globals.hh"
#include "G4VEmModel.hh"
#include "G4DataVector.hh"
#include "G4ParticleChangeForGamma.hh"

class G4ParticleDefinition;
class G4DynamicParticle;
class G4MaterialCutsCouple;
class G4Material;
class G4PhysicsFreeVector;
class G4PenelopeSamplingData;

class G4PenelopeRayleighModel : public G4VEmModel 
{

public:  
  G4PenelopeRayleighModel(const G4ParticleDefinition* p=0,
			  const G4String& processName ="PenRayleigh");
  
  virtual ~G4PenelopeRayleighModel();
  
  virtual void Initialise(const G4ParticleDefinition*, const G4DataVector&);
  virtual void InitialiseLocal(const G4ParticleDefinition*,
			       G4VEmModel *masterModel);

  virtual G4double ComputeCrossSectionPerAtom(const G4ParticleDefinition*,
                                              G4double kinEnergy,
                                              G4double Z,
                                              G4double A=0,
                                              G4double cut=0,
                                              G4double emax=DBL_MAX);

  virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
				 const G4MaterialCutsCouple*,
				 const G4DynamicParticle*,
				 G4double tmin,
				 G4double maxEnergy);
    
  void SetVerbosityLevel(G4int lev){verboseLevel = lev;};
  G4int GetVerbosityLevel(){return verboseLevel;};
  
  //Testing purposes
  void DumpFormFactorTable(const G4Material*);

protected:
  G4ParticleChangeForGamma* fParticleChange;
  const G4ParticleDefinition* fParticle;

private:
  G4PenelopeRayleighModel& operator=(const G4PenelopeRayleighModel &right);
  G4PenelopeRayleighModel(const G4PenelopeRayleighModel&);
    
  void SetParticle(const G4ParticleDefinition*);

  //Intrinsic energy limits of the model: cannot be extended by 
  //the parent process
  G4double fIntrinsicLowEnergyLimit;
  G4double fIntrinsicHighEnergyLimit;
  
  G4int verboseLevel;
  G4bool isInitialised;

  //Internal tables and manager methods
  std::map<G4int,G4PhysicsFreeVector*> *logAtomicCrossSection;
  std::map<G4int,G4PhysicsFreeVector*> *atomicFormFactor;


  G4DataVector logQSquareGrid; //log(Q^2) grid for interpolation
  std::map<const G4Material*,G4PhysicsFreeVector*> *logFormFactorTable; //log(Q^2) vs. log(F^2)
 
  G4DataVector logEnergyGridPMax; //energy grid for PMac (and originally for the x-section)
  std::map<const G4Material*,G4PhysicsFreeVector*> *pMaxTable; //E vs. Pmax

  std::map<const G4Material*,G4PenelopeSamplingData*> *samplingTable;

  //Helper methods
  void ReadDataFile(G4int); 
  void ClearTables();
  void BuildFormFactorTable(const G4Material*);
  void GetPMaxTable(const G4Material*);

  G4double GetFSquared(const G4Material*,const G4double);
  void InitializeSamplingAlgorithm(const G4Material*);

  //Used only for G4EmCalculator and Unit Tests
  G4bool fLocalTable;

};

#endif

