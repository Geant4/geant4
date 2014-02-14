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
// $Id: G4LowEWentzelVIModel.hh 74697 2013-10-19 16:15:25Z vnivanch $
//
// -------------------------------------------------------------------
//
//
// GEANT4 Class header file
//
//
// File name:     G4LowEWentzelVIModel
//
// Author:        V.Ivanchenko 
//
// Creation date: 11.02.2014 from G4WentzelVIModel
//
// Modifications:
//
// Class Description:
//
// Implementation of the model of multiple scattering for low-energy e-

// -------------------------------------------------------------------
//

#ifndef G4LowEWentzelVIModel_h
#define G4LowEWentzelVIModel_h 1

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4VMscModel.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4WentzelOKandVIxSection.hh"

class G4ParticleDefinition;
class G4LossTableManager;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4LowEWentzelVIModel : public G4VMscModel
{

public:

  G4LowEWentzelVIModel();

  virtual ~G4LowEWentzelVIModel();

  virtual void Initialise(const G4ParticleDefinition*, const G4DataVector&);

  void StartTracking(G4Track*);

  virtual G4double ComputeCrossSectionPerAtom(const G4ParticleDefinition*,
					      G4double KineticEnergy,
					      G4double AtomicNumber,
					      G4double AtomicWeight=0., 
					      G4double cut = DBL_MAX,
					      G4double emax= DBL_MAX);

  virtual G4ThreeVector& SampleScattering(const G4ThreeVector&, 
					  G4double safety);

  virtual G4double ComputeTruePathLengthLimit(const G4Track& track,
					      G4double& currentMinimalStep);

  virtual G4double ComputeGeomPathLength(G4double truePathLength);

  virtual G4double ComputeTrueStepLength(G4double geomStepLength);

  // defines low energy limit on energy transfer to atomic electron
  inline void SetFixedCut(G4double);

  // low energy limit on energy transfer to atomic electron
  inline G4double GetFixedCut() const;

  // access to cross section class
  inline G4WentzelOKandVIxSection* GetWVICrossSection();

private:

  G4double ComputeXSectionPerVolume();

  inline void SetupParticle(const G4ParticleDefinition*);

  inline void DefineMaterial(const G4MaterialCutsCouple*);

  //  hide assignment operator
  G4LowEWentzelVIModel & operator=(const  G4LowEWentzelVIModel &right);
  G4LowEWentzelVIModel(const  G4LowEWentzelVIModel&);

  G4LossTableManager*       theManager;
  G4ParticleChangeForMSC*   fParticleChange;
  G4WentzelOKandVIxSection* wokvi;

  const G4DataVector*       currentCuts;

  G4double tlimitminfix;
  G4double invsqrt12;
  G4double fixedCut;

  // cache kinematics
  G4double preKinEnergy;
  G4double tPathLength;
  G4double zPathLength;
  G4double lambdaeff;
  G4double currentRange; 

  // data for single scattering mode
  G4double xtsec;
  std::vector<G4double> xsecn;
  std::vector<G4double> prob;
  G4int    nelments;

  G4double numlimit;

  // cache material
  G4int    currentMaterialIndex;
  const G4MaterialCutsCouple* currentCouple;
  const G4Material* currentMaterial;

  // single scattering parameters
  G4double cosThetaMin;
  G4double cosThetaMax;
  G4double cosTetMaxNuc;

  // projectile
  const G4ParticleDefinition* particle;
  G4double lowEnergyLimit;

  // flags
  G4bool   inside;
  G4bool   singleScatteringMode;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline
void G4LowEWentzelVIModel::DefineMaterial(const G4MaterialCutsCouple* cup) 
{ 
  if(cup != currentCouple) {
    currentCouple = cup;
    SetCurrentCouple(cup); 
    currentMaterial = cup->GetMaterial();
    currentMaterialIndex = currentCouple->GetIndex(); 
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4LowEWentzelVIModel::SetupParticle(const G4ParticleDefinition* p)
{
  // Initialise mass and charge
  if(p != particle) {
    particle = p;
    wokvi->SetupParticle(p);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4LowEWentzelVIModel::SetFixedCut(G4double val)
{
  fixedCut = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4LowEWentzelVIModel::GetFixedCut() const
{
  return fixedCut;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4WentzelOKandVIxSection* G4LowEWentzelVIModel::GetWVICrossSection()
{
  return wokvi;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

