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
// $Id: G4PairProductionRelModel.hh 104477 2017-06-01 07:39:33Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4PairProductionRelModel
//
// Author:        Andreas Schaelicke
//
// Creation date: 02.04.2009
//
// Modifications:
//
// Class Description:
//
// Implementation of gamma convertion to e+e- in the field of a nucleus 
// relativistic approximation
// 

// -------------------------------------------------------------------
//

#ifndef G4PairProductionRelModel_h
#define G4PairProductionRelModel_h 1

#include <CLHEP/Units/PhysicalConstants.h>

#include "G4VEmModel.hh"
#include "G4PhysicsTable.hh"
#include "G4NistManager.hh"
#include "G4Log.hh"
#include "G4Exp.hh"

class G4ParticleChangeForGamma;

class G4PairProductionRelModel : public G4VEmModel
{

public:

  explicit G4PairProductionRelModel(const G4ParticleDefinition* p = nullptr, 
				    const G4String& nam = "BetheHeitlerLPM");
 
  virtual ~G4PairProductionRelModel();

  virtual void Initialise(const G4ParticleDefinition*, 
			  const G4DataVector&) override;

  virtual void InitialiseLocal(const G4ParticleDefinition*, 
			       G4VEmModel* masterModel) override;

  virtual G4double ComputeCrossSectionPerAtom(
                                const G4ParticleDefinition*,
                                      G4double kinEnergy, 
                                      G4double Z, 
                                      G4double A=0., 
                                      G4double cut=0.,
                                      G4double emax=DBL_MAX) override;

  virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
				 const G4MaterialCutsCouple*,
				 const G4DynamicParticle*,
				 G4double tmin,
				 G4double maxEnergy) override;

  virtual void SetupForMaterial(const G4ParticleDefinition*,
                                const G4Material*,G4double) override;

  // * fast inline functions *
  inline void SetCurrentElement(G4double Z);

  // set / get methods
  inline void SetLPMconstant(G4double val);
  inline G4double LPMconstant() const;

  inline void SetLPMflag(G4bool);
  inline G4bool LPMflag() const;

protected:

  // screening functions
  inline G4double Phi1(G4double delta) const;
  inline G4double Phi2(G4double delta) const;
  inline G4double ScreenFunction1(G4double ScreenVariable);
  inline G4double ScreenFunction2(G4double ScreenVariable);
  inline G4double DeltaMax() const;
  inline G4double DeltaMin(G4double) const;

  // lpm functions
  void  CalcLPMFunctions(G4double k, G4double eplus);

  G4double ComputeXSectionPerAtom(G4double totalEnergy, G4double Z);

  G4double ComputeDXSectionPerAtom(G4double eplusEnergy, G4double totalEnergy, G4double Z);
  G4double ComputeRelDXSectionPerAtom(G4double eplusEnergy, G4double totalEnergy, G4double Z);

  // hide assignment operator
  G4PairProductionRelModel & operator=
  (const G4PairProductionRelModel &right) = delete;
  G4PairProductionRelModel(const  G4PairProductionRelModel&) = delete;

  G4NistManager*              nist;

  G4ParticleDefinition*     theGamma;
  G4ParticleDefinition*     theElectron;
  G4ParticleDefinition*     thePositron;
  G4ParticleChangeForGamma* fParticleChange;

  G4double fLPMconstant;
  G4bool   fLPMflag;

  // cash
  G4double z13, z23, lnZ;
  G4double Fel, Finel, fCoulomb; 
  G4double currentZ;

  // LPM effect
  G4double lpmEnergy;
  G4double xiLPM, phiLPM, gLPM;

  // consts
  G4bool   use_completescreening;

  static const G4double xgi[8], wgi[8];
  static const G4double Fel_light[5];
  static const G4double Finel_light[5];
  static const G4double facFel;
  static const G4double facFinel;

  static const G4double preS1, logTwo, xsfactor, Egsmall, Eghigh;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline 
void G4PairProductionRelModel::SetLPMconstant(G4double val) 
{
  fLPMconstant = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline 
G4double G4PairProductionRelModel::LPMconstant() const 
{
  return fLPMconstant;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline 
void G4PairProductionRelModel::SetLPMflag(G4bool val) 
{
  fLPMflag = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline 
G4bool G4PairProductionRelModel::LPMflag() const 
{
  return fLPMflag;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void G4PairProductionRelModel::SetCurrentElement(G4double Z)
{
  if(Z != currentZ) {
    currentZ = Z;

    G4int iz = G4lrint(Z);
    z13 = nist->GetZ13(iz);
    z23 = z13*z13;
    lnZ = nist->GetLOGZ(iz);

    if (iz <= 4) {
      Fel = Fel_light[iz];  
      Finel = Finel_light[iz] ; 
    }
    else {
      Fel = facFel - lnZ/3. ;
      Finel = facFinel - 2.*lnZ/3. ;
    }
    fCoulomb=GetCurrentElement()->GetfCoulomb();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4PairProductionRelModel::Phi1(G4double delta) const
{
  return (delta > 1.) 
    ? 21.12 - 4.184*G4Log(delta+0.952)
    : 20.868 - delta*(3.242 - 0.625*delta);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4double G4PairProductionRelModel::Phi2(G4double delta) const
{
  return (delta > 1.)
    ? 21.12 - 4.184*G4Log(delta+0.952)
    : 20.209 - delta*(1.930 + 0.086*delta);
}

inline G4double G4PairProductionRelModel::ScreenFunction1(G4double ScreenVariable)
// compute the value of the screening function 3*PHI1 - PHI2
{
  return (ScreenVariable > 1.)
    ? 42.24 - 8.368*G4Log(ScreenVariable+0.952)
    : 42.392 - ScreenVariable*(7.796 - 1.961*ScreenVariable);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4double G4PairProductionRelModel::ScreenFunction2(G4double ScreenVariable)
// compute the value of the screening function 1.5*PHI1 + 0.5*PHI2
{
  return (ScreenVariable > 1.)
    ? 42.24 - 8.368*G4Log(ScreenVariable+0.952)
    : 41.405 - ScreenVariable*(5.828 - 0.8945*ScreenVariable);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4PairProductionRelModel::DeltaMax() const
{
  // k > 50 MeV
  G4double FZ = 8.*(lnZ/3. + fCoulomb);
  return G4Exp( (42.24-FZ)/8.368 ) + 0.952;
}

inline G4double G4PairProductionRelModel::DeltaMin(G4double k) const
{
  return 544.*CLHEP::electron_mass_c2/(z13*k);
}

#endif
