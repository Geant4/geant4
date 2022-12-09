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
//
// File name:     G4PAIPhotModel
//
// Author:        V. Grichine based on G4PAIModel  code for MT
//
// Creation date: 07.10.2013
//
//
// Class Description:
//
// Implementation of PAI photon model of energy loss and
// delta-electron or X-ray production by charged particles in MT framework
//
// -------------------------------------------------------------------
//

#ifndef G4PAIPhotModel_h
#define G4PAIPhotModel_h 1

#include <CLHEP/Units/PhysicalConstants.h>

#include "G4VEmModel.hh"
#include "G4VEmFluctuationModel.hh"
#include "globals.hh"
#include <vector>

class G4Region;
class G4MaterialCutsCouple;
class G4ParticleChangeForLoss;

class G4PAIPhotData;

class G4PAIPhotModel final : public G4VEmModel, public G4VEmFluctuationModel
{

public:

  explicit G4PAIPhotModel(const G4ParticleDefinition* p = nullptr, 
			  const G4String& nam = "PAI");

  ~G4PAIPhotModel() final;

  void Initialise(const G4ParticleDefinition*, const G4DataVector&) final;

  void InitialiseLocal(const G4ParticleDefinition*, 
                       G4VEmModel* masterModel) final;

  G4double MinEnergyCut(const G4ParticleDefinition*,
                        const G4MaterialCutsCouple* couple) final;

  G4double ComputeDEDXPerVolume(const G4Material*,
		                const G4ParticleDefinition*,
			        G4double kineticEnergy,
			        G4double cutEnergy) final;

  G4double CrossSectionPerVolume(const G4Material*,
			         const G4ParticleDefinition*,
				 G4double kineticEnergy,
				 G4double cutEnergy,
				 G4double maxEnergy) final;

  void SampleSecondaries(std::vector<G4DynamicParticle*>*,
				 const G4MaterialCutsCouple*,
				 const G4DynamicParticle*,
				 G4double tmin,
				 G4double maxEnergy) final;

  G4double SampleFluctuations(const G4MaterialCutsCouple*,
			      const G4DynamicParticle*,
			      const G4double, const G4double, 
                              const G4double, const G4double) final;

  G4double Dispersion(const G4Material*, const G4DynamicParticle*,
		      const G4double, const G4double, const G4double) final;

  void DefineForRegion(const G4Region* r) final;

  inline G4PAIPhotData* GetPAIPhotData();

  inline const std::vector<const G4MaterialCutsCouple*>& GetVectorOfCouples();

  inline G4double ComputeMaxEnergy(G4double scaledEnergy);

  inline void SetVerboseLevel(G4int verbose);

  // hide assignment operator 
  G4PAIPhotModel & operator=(const  G4PAIPhotModel &right) = delete;
  G4PAIPhotModel(const  G4PAIPhotModel&) = delete;

protected:

  G4double MaxSecondaryEnergy(const G4ParticleDefinition*,
				      G4double kinEnergy) final;

private:

  inline G4int FindCoupleIndex(const G4MaterialCutsCouple*);

  inline void SetParticle(const G4ParticleDefinition* p);

  G4int                       fVerbose; 

  G4PAIPhotData*             fModelData; 

  std::vector<const G4MaterialCutsCouple*> fMaterialCutsCoupleVector;
  std::vector<const G4Region*>      fPAIRegionVector;

  const G4ParticleDefinition* fParticle;
  const G4ParticleDefinition* fElectron;
  const G4ParticleDefinition* fPositron;
  G4ParticleChangeForLoss*    fParticleChange;

  G4double fMass;
  G4double fRatio;
  G4double fChargeSquare;
  G4double fLowestTcut;
};

inline G4PAIPhotData* G4PAIPhotModel::GetPAIPhotData()
{
  return fModelData;
}

inline const std::vector<const G4MaterialCutsCouple*>& 
G4PAIPhotModel::GetVectorOfCouples()
{
  return fMaterialCutsCoupleVector;
}

inline G4double G4PAIPhotModel::ComputeMaxEnergy(G4double scaledEnergy)
{
  return MaxSecondaryEnergy(fParticle, scaledEnergy/fRatio);
}

inline void G4PAIPhotModel::SetVerboseLevel(G4int verbose) 
{ 
  fVerbose=verbose; 
}

inline G4int G4PAIPhotModel::FindCoupleIndex(const G4MaterialCutsCouple* couple)
{
  G4int idx = -1;
  G4int jMatMax = (G4int)fMaterialCutsCoupleVector.size();
  for(G4int jMat = 0;jMat < jMatMax; ++jMat) { 
    if(couple == fMaterialCutsCoupleVector[jMat]) {
      idx = jMat; 
      break; 
    }
  }
  return idx;
}

inline void G4PAIPhotModel::SetParticle(const G4ParticleDefinition* p)
{
  if(fParticle != p) {
    fParticle = p;
    fMass = fParticle->GetPDGMass();
    fRatio = CLHEP::proton_mass_c2/fMass;
    G4double q = fParticle->GetPDGCharge()/CLHEP::eplus;
    fChargeSquare = q*q;
  }
}

#endif
