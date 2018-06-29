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
// GEANT4 Class header file
//
// File name:     G4AtimaEnergyLossModel
//
// Author:        Jose Luis Rodriguez Sanchez on base of ATIMA code
//
// Creation date: 16.01.2018
//
// Modifications:
//
//
// Class Description:
//
// Implementation of ATIMA model of energy loss 
// by heavy charged particles

// -------------------------------------------------------------------
//

#ifndef G4AtimaEnergyLossModel_h
#define G4AtimaEnergyLossModel_h 1

#include <CLHEP/Units/SystemOfUnits.h>

#include "G4VEmModel.hh"
#include "G4NistManager.hh"

class G4EmCorrections;
class G4ParticleChangeForLoss;

class G4AtimaEnergyLossModel : public G4VEmModel
{

public:

  explicit G4AtimaEnergyLossModel(const G4ParticleDefinition* p = nullptr,
                                  const G4String& nam = "Atima");

  virtual ~G4AtimaEnergyLossModel();

  virtual void Initialise(const G4ParticleDefinition*, 
			  const G4DataVector&) override;

  virtual G4double MinEnergyCut(const G4ParticleDefinition*,
				const G4MaterialCutsCouple* couple) override;

  virtual G4double ComputeCrossSectionPerElectron(
				 const G4ParticleDefinition*,
				 G4double kineticEnergy,
				 G4double cutEnergy,
				 G4double maxEnergy);
				 
  virtual G4double ComputeCrossSectionPerAtom(
				 const G4ParticleDefinition*,
				 G4double kineticEnergy,
				 G4double Z, G4double A,
				 G4double cutEnergy,
				 G4double maxEnergy) override;
				 				 
  virtual G4double CrossSectionPerVolume(const G4Material*,
				 const G4ParticleDefinition*,
				 G4double kineticEnergy,
				 G4double cutEnergy,
				 G4double maxEnergy) override;
				 
  virtual G4double ComputeDEDXPerVolume(const G4Material*,
					const G4ParticleDefinition*,
					G4double kineticEnergy,
					G4double) override;

  virtual G4double GetChargeSquareRatio(const G4ParticleDefinition* p,
					const G4Material* mat,
					G4double kineticEnergy) override;

  virtual G4double GetParticleCharge(const G4ParticleDefinition* p,
				     const G4Material* mat,
				     G4double kineticEnergy) override;

  virtual void CorrectionsAlongStep(const G4MaterialCutsCouple*,
				    const G4DynamicParticle*,
				    G4double&,
				    G4double&,
				    G4double) override;

  virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
				 const G4MaterialCutsCouple*,
				 const G4DynamicParticle*,
				 G4double tmin,
				 G4double maxEnergy) override;

protected:

  virtual G4double MaxSecondaryEnergy(const G4ParticleDefinition*,
				      G4double kinEnergy) override;

  inline G4double GetChargeSquareRatio() const;

  inline void SetChargeSquareRatio(G4double val);

private:

  void SetupParameters();

  inline void SetParticle(const G4ParticleDefinition* p);

  inline void SetGenericIon(const G4ParticleDefinition* p);

  G4double StoppingPower(G4double ap, G4double zp, G4double ep, G4double at, G4double zt);
  G4double Bethek_dedx_e(G4double ap,G4double zp,G4double ep,G4double at,G4double zt);
  G4double dedx_n(const G4double ap, const G4double zp, const G4double ep, const G4double at, const G4double zt);
  G4double sezi_dedx_e(const G4double zp, const G4double ep, const G4double at, const G4double zt);
  G4double sezi_p_se(const G4double energy, const G4double at, const G4double zt);
  G4double EnergyTable_interpolate(G4double xval, const G4double* y);

  // hide assignment operator
  G4AtimaEnergyLossModel & operator=(const  G4AtimaEnergyLossModel &right) = delete;
  G4AtimaEnergyLossModel(const  G4AtimaEnergyLossModel&) = delete;

  const G4ParticleDefinition* particle;
  G4ParticleDefinition*       theElectron;
  G4EmCorrections*            corr;
  G4ParticleChangeForLoss*    fParticleChange;
  G4NistManager*              nist;
  G4Pow* g4calc;

  G4double mass;
  G4double tlimit;
  G4double spin;
  G4double magMoment2;
  G4double chargeSquare;
  G4double ratio;
  G4double formfact;
  G4double corrFactor;
  G4bool   isIon;
  G4double MLN10;
  G4double atomic_mass_unit;
  G4double dedx_constant;
  G4double electron_mass;
  G4double fine_structure;
  G4double domega2dx_constant;

  static G4double stepE;
  static G4double tableE[200];
  static const G4double element_atomic_weights[110];
  static const G4double ls_coefficients_a[110][200];
  static const G4double ls_coefficients_ahi[110][200];
  static const G4double proton_stopping_coef[92][8];
  static const G4double ionisation_potentials_z[121];

  static const G4double atima_vfermi[92];
  static const G4double atima_lambda_screening[92];
  static const G4double x0[92];
  static const G4double x1[92];
  static const G4double afermi[92];
  static const G4double c[92];
  static const G4double m0[92];
  static const G4double del_0[92];

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4AtimaEnergyLossModel::SetParticle(const G4ParticleDefinition* p)
{
  if(particle != p) {
    particle = p;
    if(p->GetBaryonNumber() > 3 || p->GetPDGCharge() > CLHEP::eplus) 
      { isIon = true; } 
    SetupParameters();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4AtimaEnergyLossModel::SetGenericIon(const G4ParticleDefinition* p)
{
  if(p && p->GetParticleName() == "GenericIon") { isIon = true; }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4AtimaEnergyLossModel::GetChargeSquareRatio() const
{
  return chargeSquare;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4AtimaEnergyLossModel::SetChargeSquareRatio(G4double val)
{
  chargeSquare = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
