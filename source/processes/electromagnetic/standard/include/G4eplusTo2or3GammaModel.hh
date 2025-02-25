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
// File name:     G4eplusTo2or3GammaModel
//
// Author:        Vladimir Ivanchenko and Omrane Kadri 
//
// Creation date: 29.03.2018
//
//
// Class Description:
//
// Implementation of e+ annihilation into 2 or 3 gamma on fly
// Annihilation at rest is sampled by G4VPositronAtRestModel
// Implementation of 3-gamma annihilation is performed by
// G4eplusTo3GammaOKVIModel. Cross section of both models
// depend on cut parameter fDelta, which defines relative low-energy
// limit on 3d gamma energy. For computation of the cross section
// next to leading order radiative corrections are taken into account,
// atomic effects at low-energy are not considered.
//
// V.N.Baier, V.S. Fadin, V.A. Khose, E.A. Kuraev,
// Physics Reports 78 (1981) 293-393.
//
// -------------------------------------------------------------------
//

#ifndef G4eplusTo2or3GammaModel_h
#define G4eplusTo2or3GammaModel_h 1

#include "G4VEmModel.hh"

class G4eplusTo3GammaOKVIModel;
class G4ParticleChangeForGamma;
class G4PhysicsVector;
class G4DataVector;

class G4eplusTo2or3GammaModel : public G4VEmModel
{

public:

  G4eplusTo2or3GammaModel();

  ~G4eplusTo2or3GammaModel() override;

  void Initialise(const G4ParticleDefinition*, const G4DataVector&) override;

  G4double ComputeCrossSectionPerElectron(G4double kinEnergy); 
  
  G4double ComputeCrossSectionPerAtom(
                                 const G4ParticleDefinition*,
                                 G4double kinEnergy, 
                                 G4double Z, 
                                 G4double A = 0., 
                                 G4double cutEnergy = 0.,
                                 G4double maxEnergy = DBL_MAX) override;

  G4double CrossSectionPerVolume(const G4Material*,
				 const G4ParticleDefinition*,
				 G4double kineticEnergy,
				 G4double cutEnergy = 0.0,
				 G4double maxEnergy = DBL_MAX) override;

  void SampleSecondaries(std::vector<G4DynamicParticle*>*,
				 const G4MaterialCutsCouple*,
				 const G4DynamicParticle*,
				 G4double tmin = 0.0,
				 G4double maxEnergy = DBL_MAX) override;

  void SetDelta(G4double val) { if(val > 0.0) { fDeltaMin = val; } };

  // hide assignment operator
  G4eplusTo2or3GammaModel & operator=
  (const  G4eplusTo2or3GammaModel &right) = delete;
  G4eplusTo2or3GammaModel(const  G4eplusTo2or3GammaModel&) = delete;

private:

  const G4ParticleDefinition* theGamma;
  G4ParticleChangeForGamma* fParticleChange;
  G4eplusTo3GammaOKVIModel* f3GModel;

  G4double fDeltaMin; // fixed minimal relative limit
  G4double fDelta;    // running limit - function of energy
  G4double fGammaTh;  // 3-gamma annihilation low-energy limit

  static G4PhysicsVector* fCrossSection;
  static G4PhysicsVector* f3GProbability;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
