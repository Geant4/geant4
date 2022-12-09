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
// Author: Luciano Pandola
//
// History:
// -----------
// 08 Jan 2010   L. Pandola   1st implementation. 
// 25 May 2011   L. Pandola   Renamed (make v2008 as default Penelope)
// 10 Jun 2011   L. Pandola   Migrated to the new AtomDeexcitation interface
// 18 Sep 2013   L. Pandola   Migrated to the MT paradigm
//
// -------------------------------------------------------------------
//
// Class description:
// Low Energy Electromagnetic Physics, Photo-electric effect
// with Penelope Model, version 2008
// -------------------------------------------------------------------

#ifndef G4PENELOPEPHOTOELECTRICMODEL_HH
#define G4PENELOPEPHOTOELECTRICMODEL_HH 1

#include "globals.hh"
#include "G4VEmModel.hh"
#include "G4DataVector.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4AtomicTransitionManager.hh"
#include "G4VAtomDeexcitation.hh"

class G4ParticleDefinition;
class G4DynamicParticle;
class G4MaterialCutsCouple;
class G4Material;

class G4PenelopePhotoElectricModel : public G4VEmModel 
{
public:
  explicit G4PenelopePhotoElectricModel(const G4ParticleDefinition* p=nullptr,
					const G4String& processName ="PenPhotoElec");
  virtual ~G4PenelopePhotoElectricModel();

  void Initialise(const G4ParticleDefinition*, const G4DataVector&) override;
  void InitialiseLocal(const G4ParticleDefinition*,
                               G4VEmModel *masterModel) override;

  G4double ComputeCrossSectionPerAtom(
				      const G4ParticleDefinition*,
				      G4double kinEnergy,
				      G4double Z,
				      G4double A=0,
				      G4double cut=0,
				      G4double emax=DBL_MAX) override;

  void SampleSecondaries(std::vector<G4DynamicParticle*>*,
			 const G4MaterialCutsCouple*,
			 const G4DynamicParticle*,
			 G4double tmin,
			 G4double maxEnergy) override;

  void SetVerbosityLevel(G4int lev){fVerboseLevel = lev;};
  G4int GetVerbosityLevel(){return fVerboseLevel;};

  //testing purposes
  std::size_t GetNumberOfShellXS(G4int);
  G4double GetShellCrossSection(G4int Z,std::size_t shellID,G4double energy);

  G4PenelopePhotoElectricModel & operator=(const G4PenelopePhotoElectricModel &right) = delete;
  G4PenelopePhotoElectricModel(const G4PenelopePhotoElectricModel&) = delete;

protected:
  G4ParticleChangeForGamma* fParticleChange;
  const G4ParticleDefinition* fParticle;

private:
  void SetParticle(const G4ParticleDefinition*);
  G4double SampleElectronDirection(G4double energy);
  void ReadDataFile(G4int Z);
  std::size_t SelectRandomShell(G4int Z,G4double energy);
  G4String WriteTargetShell(std::size_t shellID);

  //For each Z, the PhysicsTable contains nShell+1 physics vectors
  //with log(E) vs. log(XS)
  //Element [0] of the table is the total XS, element [iS] is the 
  //partial cross section for shell iS-1
  static const G4int fMaxZ =99;
  static G4PhysicsTable* fLogAtomicShellXS[fMaxZ+1];

  G4VAtomDeexcitation*             fAtomDeexcitation;
  const G4AtomicTransitionManager* fTransitionManager;

  //Intrinsic energy limits of the model: cannot be extended by the parent process
  G4double fIntrinsicLowEnergyLimit;
  G4double fIntrinsicHighEnergyLimit;
  G4int fVerboseLevel;
  G4bool fIsInitialised; 
  //Used only for G4EmCalculator and Unit Tests
  G4bool fLocalTable;
};

#endif

