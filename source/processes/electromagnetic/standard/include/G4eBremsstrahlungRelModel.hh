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
// $Id: G4eBremsstrahlungRelModel.hh 104373 2017-05-29 09:56:39Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4eBremsstrahlungRelModel
//                extention of standard G4eBremsstrahlungModel
//
// Author:        Andreas Schaelicke 
//
// Creation date: 28.03.2008
//
// Modifications:
//
//
// Class Description:
//
// Implementation of energy loss for gamma emission by electrons and
// positrons including an improved version of the LPM effect

// -------------------------------------------------------------------
//

#ifndef G4eBremsstrahlungRelModel_h
#define G4eBremsstrahlungRelModel_h 1

#include "G4VEmModel.hh"
#include "G4NistManager.hh"
#include "G4Exp.hh"
#include "G4Log.hh"

class G4ParticleChangeForLoss;
class G4PhysicsVector;

class G4eBremsstrahlungRelModel : public G4VEmModel
{

public:

  explicit G4eBremsstrahlungRelModel(const G4ParticleDefinition* p = 0, 
				     const G4String& nam = "eBremLPM");

  virtual ~G4eBremsstrahlungRelModel();

  virtual void Initialise(const G4ParticleDefinition*, const G4DataVector&) override;

  virtual void InitialiseLocal(const G4ParticleDefinition*,
			       G4VEmModel* masterModel) override;

  virtual G4double ComputeDEDXPerVolume(const G4Material*,
					const G4ParticleDefinition*,
					G4double kineticEnergy,
					G4double cutEnergy) override;
					
  virtual G4double ComputeCrossSectionPerAtom(const G4ParticleDefinition*,
					      G4double tkin, 
					      G4double Z,   G4double,
					      G4double cutEnergy,
					      G4double maxEnergy = DBL_MAX) override;
  
  virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
				 const G4MaterialCutsCouple*,
				 const G4DynamicParticle*,
				 G4double cutEnergy,
				 G4double maxEnergy) override;

  virtual void SetupForMaterial(const G4ParticleDefinition*,
                                const G4Material*,G4double) override;

  virtual G4double MinPrimaryEnergy(const G4Material*,
				    const G4ParticleDefinition*,
				    G4double cut) override;

  inline void SetLPMconstant(G4double val);
  inline G4double LPMconstant() const;

  inline void SetLowestKinEnergy(G4double);
  inline G4double LowestKinEnergy() const;


protected:

  virtual G4double ComputeDXSectionPerAtom(G4double gammaEnergy);

  // * fast inline functions *
  inline void SetCurrentElement(G4int);

private:

  void InitialiseConstants();

  void CalcLPMFunctions(G4double gammaEnergy);

  G4double ComputeBremLoss(G4double cutEnergy);

  G4double ComputeXSectionPerAtom(G4double cutEnergy);

  G4double ComputeRelDXSectionPerAtom(G4double gammaEnergy);

  void SetParticle(const G4ParticleDefinition* p);

  inline G4double Phi1(G4double,G4double);
  inline G4double Phi1M2(G4double,G4double);
  inline G4double Psi1(G4double,G4double);
  inline G4double Psi1M2(G4double,G4double);

  // hide assignment operator
  G4eBremsstrahlungRelModel & operator=(const  G4eBremsstrahlungRelModel &right) = delete;
  G4eBremsstrahlungRelModel(const  G4eBremsstrahlungRelModel&) = delete;

protected:

  G4NistManager*              nist;
  const G4ParticleDefinition* particle;
  G4ParticleDefinition*       theGamma;
  G4ParticleChangeForLoss*    fParticleChange;

  G4double bremFactor;

  // cash
  G4double particleMass;
  G4double kinEnergy;
  G4double totalEnergy;
  G4double densityFactor;
  G4double densityCorr;

  // scattering off electrons
  G4double nucTerm;
  G4double sumTerm;

  G4bool scatOffElectron;
  G4bool isElectron;

  G4int  currentZ;

private:

  static const G4double xgi[8], wgi[8];
  static const G4double Fel_light[5];
  static const G4double Finel_light[5];

  // consts
  G4double lowestKinEnergy;
  G4double fMigdalConstant;
  G4double fLPMconstant;
  G4double energyThresholdLPM;
  G4double facFel, facFinel;
  G4double preS1,logTwo;

  // cash
  G4double z13, z23, lnZ;
  G4double Fel, Finel, fCoulomb, fMax; 

  // LPM effect
  G4double lpmEnergy;
  G4double xiLPM, phiLPM, gLPM;

  // critical gamma energies
  G4double klpm, kp;

  // flags
  G4bool   use_completescreening;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void G4eBremsstrahlungRelModel::SetCurrentElement(G4int Z)
{
  if(Z != currentZ) {
    currentZ = Z;

    z13 = nist->GetZ13(Z);
    z23 = z13*z13;
    lnZ = nist->GetLOGZ(Z);

    if (Z <= 4) {
      Fel = Fel_light[Z];  
      Finel = Finel_light[Z] ; 
    }
    else {
      G4double lnzt = lnZ/3.;
      Fel = facFel - lnzt;
      Finel = facFinel - 2*lnzt;
    }

    fCoulomb = GetCurrentElement()->GetfCoulomb();
    G4double xz = 1.0/(G4double)Z;
    fMax = Fel-fCoulomb + Finel*xz  + (1. + xz)/12.;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


inline G4double G4eBremsstrahlungRelModel::Phi1(G4double gg, G4double)
{
  //       Thomas-Fermi FF from Tsai, eq.(3.38) for Z>=5
  return 20.863 - 2.*G4Log(1. + sqr(0.55846*gg) )
    - 4.*( 1. - 0.6*G4Exp(-0.9*gg) - 0.4*G4Exp(-1.5*gg) );
}

inline G4double G4eBremsstrahlungRelModel::Phi1M2(G4double gg, G4double)
{
  //       Thomas-Fermi FF from Tsai, eq. (3.39) for Z>=5
  // return Phi1(gg,Z) - 
  return 2./(3.*(1. + 6.5*gg +6.*gg*gg) );
}

inline G4double G4eBremsstrahlungRelModel::Psi1(G4double eps, G4double)
{
  //       Thomas-Fermi FF from Tsai, eq.(3.40) for Z>=5 
  return 28.340 - 2.*G4Log(1. + sqr(3.621*eps) )
    - 4.*( 1. - 0.7*G4Exp(-8*eps) - 0.3*G4Exp(-29.*eps) );
}

inline G4double G4eBremsstrahlungRelModel::Psi1M2(G4double eps, G4double)
{
  //       Thomas-Fermi FF from Tsai, eq. (3.41) for Z>=5
  return  2./(3.*(1. + 40.*eps +400.*eps*eps) );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline 
void G4eBremsstrahlungRelModel::SetLPMconstant(G4double val) 
{
  fLPMconstant = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline 
G4double G4eBremsstrahlungRelModel::LPMconstant() const 
{
  return fLPMconstant;
}

inline void G4eBremsstrahlungRelModel::SetLowestKinEnergy(G4double val)
{
  lowestKinEnergy = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4eBremsstrahlungRelModel::LowestKinEnergy() const
{
  return lowestKinEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


#endif
