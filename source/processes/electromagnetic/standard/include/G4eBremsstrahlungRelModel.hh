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
// $Id: G4eBremsstrahlungRelModel.hh,v 1.2 2008-08-13 16:08:50 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

class G4ParticleChangeForLoss;

class G4eBremsstrahlungRelModel : public G4VEmModel
{

public:

  G4eBremsstrahlungRelModel(const G4ParticleDefinition* p = 0, 
			    const G4String& nam = "eBremRel");

  virtual ~G4eBremsstrahlungRelModel();

  virtual void Initialise(const G4ParticleDefinition*, const G4DataVector&);

  G4double MinEnergyCut(const G4ParticleDefinition*, 
			const G4MaterialCutsCouple*);

  virtual G4double ComputeDEDXPerVolume(const G4Material*,
					const G4ParticleDefinition*,
					G4double kineticEnergy,
					G4double cutEnergy);
					
  virtual G4double ComputeCrossSectionPerAtom(const G4ParticleDefinition*,
					      G4double tkin, 
					      G4double Z,   G4double,
					      G4double cutEnergy,
					      G4double maxEnergy);
  
  virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
				 const G4MaterialCutsCouple*,
				 const G4DynamicParticle*,
				 G4double cutEnergy,
				 G4double maxEnergy);

  virtual void SetupForMaterial(const G4ParticleDefinition*,
                                const G4Material*);

  inline void SetEnergyThreshold(G4double val);
  inline G4double EnergyThreshold() const;

  inline void SetHydrogenEnergyThreshold(G4double val);
  inline G4double HydrogenEnergyThreshold() const;

protected:

  inline G4double MaxSecondaryEnergy(const G4ParticleDefinition*,
				     G4double kineticEnergy);

private:

  G4double ComputeBremLoss(G4double cutEnergy);

  G4double ComputeXSectionPerAtom(G4double cutEnergy);

  G4double ComputeDXSectionPerAtom(G4double gammaEnergy);

  G4double ComputeRelDXSectionPerAtom(G4double gammaEnergy);

  void SetParticle(const G4ParticleDefinition* p);

  inline void SetCurrentElement(G4double Z);

  // hide assignment operator
  G4eBremsstrahlungRelModel & operator=(const  G4eBremsstrahlungRelModel &right);
  G4eBremsstrahlungRelModel(const  G4eBremsstrahlungRelModel&);

protected:

  G4NistManager*              nist;
  const G4ParticleDefinition* particle;
  G4ParticleDefinition*       theGamma;
  G4ParticleChangeForLoss*    fParticleChange;

  static G4double xgi[8], wgi[8];

  G4double minThreshold;

  // cash
  G4double particleMass;
  G4double kinEnergy;
  G4double totalEnergy;
  G4double currentZ;
  G4double z13;
  G4double z23;
  G4double lnZ;
  G4double densityFactor;
  G4double densityCorr;
  G4double lpmEnergy;

  G4bool   isElectron;

private:

  G4double highKinEnergy;
  G4double lowKinEnergy;
  G4double MigdalConstant;
  G4double LPMconstant;
  G4double bremFactor;
  G4double highEnergyTh;
  G4double hydrogenEnergyTh;

  G4bool   isInitialised;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void G4eBremsstrahlungRelModel::SetCurrentElement(G4double Z)
{
  if(Z != currentZ) {
    currentZ = Z;
    G4int iz = G4int(Z);
    z13 = nist->GetZ13(iz);
    z23 = z13*z13;
    lnZ = nist->GetLOGZ(iz);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline 
G4double G4eBremsstrahlungRelModel::MaxSecondaryEnergy(
                                 const G4ParticleDefinition*,
    				       G4double kineticEnergy)
{
  return kineticEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline 
void G4eBremsstrahlungRelModel::SetEnergyThreshold(G4double val) 
{
  highEnergyTh = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline 
G4double G4eBremsstrahlungRelModel::EnergyThreshold() const 
{
  return highEnergyTh;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline 
void G4eBremsstrahlungRelModel::SetHydrogenEnergyThreshold(G4double val) 
{
  hydrogenEnergyTh = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline 
G4double G4eBremsstrahlungRelModel::HydrogenEnergyThreshold() const 
{
  return hydrogenEnergyTh;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
