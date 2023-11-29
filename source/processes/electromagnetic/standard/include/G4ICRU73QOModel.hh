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
// File name:     G4ICRU73QOModel
//
// Author:        Alexander Bagulya
//
// Creation date: 21.05.2010
//
// Modifications:
//
//
// Class Description:
//
// Quantum Harmonic Oscillator Model for energy loss using atomic shell 
// structure of atoms taking into account Q^2 (main for projectile charge Q), 
// Q^3 and Q^4 terms for computation of energy loss due to binary collisions. 
// Can be applied on heavy negatively charged particles for the energy interval 
// 10 keV - 10 MeV scaled to the proton mass.
//
// Used data and formula of 
// 1. G4QAOLowEnergyLoss class, S.Chauvie, P.Nieminen, M.G.Pia. IEEE Trans. 
//    Nucl. Sci. 54 (2007) 578.
// 2. ShellStrength and ShellEnergy from ICRU'73 Report 2005,
// 3. Data for Ta (Z=73) from P.Sigmund, A.Shinner. Eur. Phys. J. D15 (2001) 
//    165-172
//
// -------------------------------------------------------------------
//

#ifndef G4ICRU73QOModel_h
#define G4ICRU73QOModel_h 1

#include <CLHEP/Units/PhysicalConstants.h>

#include "G4VEmModel.hh"
#include "G4AtomicShells.hh"
#include "G4DensityEffectData.hh"

class G4ParticleChangeForLoss;

class G4ICRU73QOModel : public G4VEmModel
{

public:

  explicit G4ICRU73QOModel(const G4ParticleDefinition* p = nullptr,
                           const G4String& nam = "ICRU73QO");

  ~G4ICRU73QOModel() = default;

  void Initialise(const G4ParticleDefinition*, const G4DataVector&) override;

  G4double ComputeCrossSectionPerElectron(
                                 const G4ParticleDefinition*,
                                 G4double kineticEnergy,
                                 G4double cutEnergy,
                                 G4double maxEnergy);
                                 
  G4double ComputeCrossSectionPerAtom(
                                 const G4ParticleDefinition*,
                                 G4double kineticEnergy,
                                 G4double Z, G4double A,
                                 G4double cutEnergy,
                                 G4double maxEnergy) override;
                                                                  
  G4double CrossSectionPerVolume(const G4Material*,
                                 const G4ParticleDefinition*,
                                 G4double kineticEnergy,
                                 G4double cutEnergy,
                                 G4double maxEnergy) override;
                                 
  G4double ComputeDEDXPerVolume(const G4Material*,
				const G4ParticleDefinition*,
				G4double kineticEnergy,
				G4double) override;

  void SampleSecondaries(std::vector<G4DynamicParticle*>*,
			 const G4MaterialCutsCouple*,
			 const G4DynamicParticle*,
			 G4double tmin,
			 G4double maxEnergy) override;

  // add correction to energy loss and compute non-ionizing energy loss
  void CorrectionsAlongStep(const G4MaterialCutsCouple*,
			    const G4DynamicParticle*,
			    const G4double& length,
			    G4double& eloss) override;

  // hide assignment operator
  G4ICRU73QOModel & operator=(const  G4ICRU73QOModel &right) = delete;
  G4ICRU73QOModel(const  G4ICRU73QOModel&) = delete;
 
protected:

  G4double MaxSecondaryEnergy(const G4ParticleDefinition*,
			      G4double kinEnergy) final;

private:

  inline void SetParticle(const G4ParticleDefinition* p);
  inline void SetLowestKinEnergy(G4double val);

  G4double DEDX(const G4Material* material, G4double kineticEnergy);

  G4double DEDXPerElement(G4int Z, G4double kineticEnergy);

  // get number of shell, energy and oscillator strengths for material
  G4int GetNumberOfShells(G4int Z) const;

  G4double GetShellEnergy(G4int Z, G4int nbOfTheShell) const; 
  G4double GetOscillatorEnergy(G4int Z, G4int nbOfTheShell) const; 
  G4double GetShellStrength(G4int Z, G4int nbOfTheShell) const;

  // calculate stopping number for L's term
  G4double GetL0(G4double normEnergy) const;
  // terms in Z^2
  G4double GetL1(G4double normEnergy) const;
  // terms in Z^3
  G4double GetL2(G4double normEnergy) const;
  // terms in Z^4
  
  const G4ParticleDefinition* particle;
  G4ParticleDefinition*       theElectron;   
  G4ParticleChangeForLoss*    fParticleChange;
  G4DensityEffectData*        denEffData;

  G4double mass;
  G4double charge;
  G4double chargeSquare;
  G4double massRate;
  G4double ratio;
  G4double lowestKinEnergy;

  G4bool   isInitialised;

  // Z of element at now avaliable for the model
  static const G4int NQOELEM  = 26;
  static const G4int NQODATA  = 130;
  static const G4int ZElementAvailable[NQOELEM];
  
  // number, energy and oscillator strengths
  // for an harmonic oscillator model of material
  static const G4int startElemIndex[NQOELEM];
  static const G4int nbofShellsForElement[NQOELEM];
  static const G4double ShellEnergy[NQODATA];
  static const G4double SubShellOccupation[NQODATA];  // Z * ShellStrength

  G4int indexZ[100];

  //  variable for calculation of stopping number of L's term
  static const G4double L0[67][2];
  static const G4double L1[22][2];
  static const G4double L2[14][2];
  
  G4int sizeL0;
  G4int sizeL1;
  G4int sizeL2;

  static const G4double factorBethe[99];
  
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void G4ICRU73QOModel::SetParticle(const G4ParticleDefinition* p)
{
  particle = p;
  mass = particle->GetPDGMass();
  charge = particle->GetPDGCharge()/CLHEP::eplus;
  chargeSquare = charge*charge;
  massRate     = mass/CLHEP::proton_mass_c2;
  ratio = CLHEP::electron_mass_c2/mass;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void G4ICRU73QOModel::SetLowestKinEnergy(G4double val)
{
  lowestKinEnergy = val;
}

#endif
