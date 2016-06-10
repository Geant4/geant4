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
// $Id: G4EmCalculator.hh 83007 2014-07-24 14:46:57Z gcosmo $
//
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4EmCalculator
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 27.06.2004
//
// Modifications:
// 17.11.2004 Change signature of methods, add new methods (V.Ivanchenko)
// 11.01.2006 Add GetCSDARange (V.Ivanchenko)
// 26.01.2006 Rename GetRange -> GetRangeFromRestricteDEDX (V.Ivanchenko)
// 22.03.2006 Add ComputeElectronicDEDX and ComputeTotalDEDX (V.Ivanchenko)
// 29.09.2006 Add member loweModel (V.Ivanchenko)
// 15.03.2007 Add ComputeEnergyCutFromRangeCut methods (V.Ivanchenko)
//
// Class Description:
//
// Provide access to dE/dx and cross sections

// -------------------------------------------------------------------
//

#ifndef G4EmCalculator_h
#define G4EmCalculator_h 1

#include <vector>
#include "globals.hh"
#include "G4DataVector.hh"
#include "G4DynamicParticle.hh"
#include "G4VAtomDeexcitation.hh"

class G4LossTableManager;
class G4Material;
class G4MaterialCutsCouple;
class G4ParticleDefinition;
class G4PhysicsTable;
class G4VEmModel;
class G4VEnergyLossProcess;
class G4VEmProcess;
class G4VMultipleScattering;
class G4VProcess;
class G4ionEffectiveCharge;
class G4Region;
class G4Element;
class G4EmCorrections;
class G4IonTable;

class G4EmCalculator
{

public:

  G4EmCalculator();

  ~G4EmCalculator();

  //===========================================================================
  // Methods to access precalculated dE/dx and cross sections
  // Materials should exist in the list of the G4MaterialCutsCouple
  //===========================================================================

  G4double GetDEDX(G4double kinEnergy, const G4ParticleDefinition*, 
		   const G4Material*,
                   const G4Region* r = 0);
  G4double GetDEDX(G4double kinEnergy, const G4String& part, 
		   const G4String& mat,
                   const G4String& s = "world");

  G4double GetRangeFromRestricteDEDX(G4double kinEnergy, 
				     const G4ParticleDefinition*, 
				     const G4Material*,
				     const G4Region* r = 0);
  G4double GetRangeFromRestricteDEDX(G4double kinEnergy, const G4String& part, 
				     const G4String& mat,
				     const G4String& s = "world");

  G4double GetCSDARange(G4double kinEnergy, const G4ParticleDefinition*, 
			const G4Material*,
			const G4Region* r = 0);
  G4double GetCSDARange(G4double kinEnergy, const G4String& part, 
			const G4String& mat,
			const G4String& s = "world");

  G4double GetRange(G4double kinEnergy, const G4ParticleDefinition*, 
			const G4Material*,
			const G4Region* r = 0);
  G4double GetRange(G4double kinEnergy, const G4String& part, 
			const G4String& mat,
			const G4String& s = "world");

  G4double GetKinEnergy(G4double range, const G4ParticleDefinition*, 
			const G4Material*,
			const G4Region* r = 0);
  G4double GetKinEnergy(G4double range, const G4String& part, 
			const G4String& mat,
			const G4String& s = "world");

  G4double GetCrossSectionPerVolume(
                   G4double kinEnergy, const G4ParticleDefinition*,
                   const G4String& processName,  const G4Material*,
		   const G4Region* r = 0);
  G4double GetCrossSectionPerVolume(
                   G4double kinEnergy, const G4String& part, const G4String& proc,
                   const G4String& mat, const G4String& s = "world");

  G4double GetShellIonisationCrossSectionPerAtom(
                   const G4String& part, G4int Z, 
		   G4AtomicShellEnumerator shell,
                   G4double kinEnergy);

  G4double GetMeanFreePath(G4double kinEnergy, const G4ParticleDefinition*,
                   const G4String& processName,  const G4Material*,
		   const G4Region* r = 0);
  G4double GetMeanFreePath(G4double kinEnergy, const G4String& part, 
			   const G4String& proc,
                   const G4String& mat, const G4String& s = "world");

  void PrintDEDXTable(const G4ParticleDefinition*);

  void PrintRangeTable(const G4ParticleDefinition*);

  void PrintInverseRangeTable(const G4ParticleDefinition*);

  //===========================================================================
  // Methods to calculate dE/dx and cross sections "on fly"
  // Existing tables and G4MaterialCutsCouples are not used
  //===========================================================================

  G4double ComputeDEDX(G4double kinEnergy, const G4ParticleDefinition*,
                       const G4String& processName,  const G4Material*,
		       G4double cut = DBL_MAX);
  G4double ComputeDEDX(G4double kinEnergy, const G4String& part, 
		       const G4String& proc,
                       const G4String& mat, G4double cut = DBL_MAX);

  G4double ComputeElectronicDEDX(G4double kinEnergy, 
				 const G4ParticleDefinition*,
				 const G4Material* mat, G4double cut = DBL_MAX);
  G4double ComputeElectronicDEDX(G4double kinEnergy, const G4String& part,
				 const G4String& mat, G4double cut = DBL_MAX);

  G4double ComputeNuclearDEDX(G4double kinEnergy, const G4ParticleDefinition*, 
			      const G4Material*);
  G4double ComputeNuclearDEDX(G4double kinEnergy, const G4String& part, 
			      const G4String& mat);

  G4double ComputeTotalDEDX(G4double kinEnergy, const G4ParticleDefinition*, 
			    const G4Material*, G4double cut = DBL_MAX);
  G4double ComputeTotalDEDX(G4double kinEnergy, const G4String& part, 
			    const G4String& mat, G4double cut = DBL_MAX);

  G4double ComputeCrossSectionPerVolume(
                       G4double kinEnergy, const G4ParticleDefinition*,
                       const G4String& processName,  const G4Material*,
		       G4double cut = 0.0);
  G4double ComputeCrossSectionPerVolume(
                       G4double kinEnergy, const G4String& part, 
		       const G4String& proc,
                       const G4String& mat, G4double cut = 0.0);

  G4double ComputeCrossSectionPerAtom(
                       G4double kinEnergy, const G4ParticleDefinition*,
                       const G4String& processName, G4double Z, G4double A,
		       G4double cut = 0.0);
  G4double ComputeCrossSectionPerAtom(
                       G4double kinEnergy, const G4String& part,
                       const G4String& processName, const G4Element*,
		       G4double cut = 0.0);

  G4double ComputeGammaAttenuationLength(G4double kinEnergy, 
					 const G4Material*);

  G4double ComputeShellIonisationCrossSectionPerAtom(
                   const G4String& part, G4int Z, 
		   G4AtomicShellEnumerator shell,
                   G4double kinEnergy,
                   const G4Material* mat = 0);

  G4double ComputeMeanFreePath(
                       G4double kinEnergy, const G4ParticleDefinition*,
                       const G4String& processName,  const G4Material*,
		       G4double cut = 0.0);
  G4double ComputeMeanFreePath(
                       G4double kinEnergy, const G4String&, const G4String&,
                       const G4String& processName, G4double cut = 0.0);

  G4double ComputeEnergyCutFromRangeCut(
                       G4double range, const G4ParticleDefinition*,
		       const G4Material*);
  G4double ComputeEnergyCutFromRangeCut(
                       G4double range, const G4String&,
		       const G4String&);

  //===========================================================================
  // Methods to access particles, materials, regions, processes
  //===========================================================================

  const G4ParticleDefinition* FindParticle(const G4String&);

  const G4ParticleDefinition* FindIon(G4int Z, G4int A);

  const G4Material* FindMaterial(const G4String&);

  const G4Region* FindRegion(const G4String&);

  const G4MaterialCutsCouple* FindCouple(const G4Material*, 
					 const G4Region* r = 0);

  G4VProcess* FindProcess(const G4ParticleDefinition* part,
			  const G4String& processName);

  void SetVerbose(G4int val);

  //===========================================================================
  // Private methods 
  //===========================================================================

private:

  G4bool UpdateParticle(const G4ParticleDefinition*, G4double kinEnergy);

  G4bool UpdateCouple(const G4Material*, G4double cut);

  void FindLambdaTable(const G4ParticleDefinition*, 
		       const G4String& processName,
		       G4double kinEnergy);

  G4bool FindEmModel(const G4ParticleDefinition*, 
                     const G4String& processName,
                           G4double kinEnergy);

  G4VEnergyLossProcess* FindEnergyLossProcess(const G4ParticleDefinition*);

  G4VEnergyLossProcess* FindEnLossProcess(const G4ParticleDefinition*,
					  const G4String& processName);

  G4VEmProcess* FindDiscreteProcess(const G4ParticleDefinition*,
				    const G4String& processName);

  G4VMultipleScattering* FindMscProcess(const G4ParticleDefinition*,
					const G4String& processName);

  G4bool ActiveForParticle(const G4ParticleDefinition* part,
			   G4VProcess* proc);

  G4EmCalculator & operator=(const  G4EmCalculator &right);
  G4EmCalculator(const  G4EmCalculator&);

  std::vector<const G4Material*>            localMaterials;
  std::vector<const G4MaterialCutsCouple*>  localCouples;

  G4LossTableManager*          manager;
  G4IonTable*                  ionTable;
  G4EmCorrections*             corr; 
  G4DataVector                 localCuts;
  G4int                        nLocalMaterials;

  G4int                        verbose;

  // cash
  G4int                        currentCoupleIndex;
  const G4MaterialCutsCouple*  currentCouple;
  const G4Material*            currentMaterial;
  const G4ParticleDefinition*  currentParticle;
  const G4ParticleDefinition*  lambdaParticle;
  const G4ParticleDefinition*  baseParticle;
  const G4PhysicsTable*        currentLambda;
        G4VEmModel*            currentModel;
        G4VEmModel*            loweModel;
        G4VEnergyLossProcess*  currentProcess;

  const G4ParticleDefinition*  theGenericIon;
  G4ionEffectiveCharge*        ionEffCharge;
  G4DynamicParticle            dynParticle;

  G4String                     currentName;
  G4String                     lambdaName;
  G4double                     currentCut;
  G4double                     chargeSquare;
  G4double                     massRatio;
  G4double                     mass;
  G4bool                       isIon;
  G4bool                       isApplicable;

  G4String                     currentParticleName;
  G4String                     currentMaterialName;
  G4String                     currentProcessName;
};

//....oooOO0OOooo.......oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
