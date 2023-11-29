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
// File name:     G4MuBremsstrahlungModel
//
// Author:        Vladimir Ivanchenko on base of Laszlo Urban code
// 
// Creation date: 18.05.2002
//
// Modifications:
//
// 23-12-02 Change interface in order to move to cut per region (V.Ivanchenko)
// 27-01-03 Make models region aware (V.Ivanchenko)
// 13-02-03 Add name (V.Ivanchenko)
// 10-02-04 Add lowestKinEnergy (V.Ivanchenko)
// 08-04-05 Major optimisation of internal interfaces (V.Ivantchenko)
// 13-02-06 Add ComputeCrossSectionPerAtom (mma)
// 11-10-07 Add ignoreCut flag (V.Ivanchenko) 
// 28-02-08 Reorganized protected methods and members (V.Ivanchenko) 
// 06-03-08 Remove obsolete methods and members (V.Ivanchenko) 
// 31-05-13 Use element selectors instead of local data (V.Ivanchenko)
//
// Class Description:
//
// Implementation of bremssrahlung by muons
// A.G. Bogdanov et al., IEEE Trans. Nuc. Sci., Vol.53, No.2, 2006
//
// -------------------------------------------------------------------
//

#ifndef G4MuBremsstrahlungModel_h
#define G4MuBremsstrahlungModel_h 1

#include "G4VEmModel.hh"
#include "G4NistManager.hh"
#include "G4Exp.hh"

class G4Element;
class G4ParticleChangeForLoss;

class G4MuBremsstrahlungModel : public G4VEmModel
{

public:

  explicit G4MuBremsstrahlungModel(const G4ParticleDefinition* p = nullptr,
                                   const G4String& nam = "MuBrem");

  ~G4MuBremsstrahlungModel() = default;

  void Initialise(const G4ParticleDefinition*, const G4DataVector&) override;

  void InitialiseLocal(const G4ParticleDefinition*, 
                       G4VEmModel* masterModel) override;

  G4double MinEnergyCut(const G4ParticleDefinition*,
                        const G4MaterialCutsCouple*) override;
			      
  G4double ComputeCrossSectionPerAtom(
				 const G4ParticleDefinition*,
				 G4double kineticEnergy,
				 G4double Z, G4double A,
				 G4double cutEnergy,
				 G4double maxEnergy) override;
				 			       
  G4double ComputeDEDXPerVolume(const G4Material*,
                                const G4ParticleDefinition*,
                                G4double kineticEnergy,
                                G4double cutEnergy) override;
			      
  void SampleSecondaries(std::vector<G4DynamicParticle*>*,
                         const G4MaterialCutsCouple*,
			 const G4DynamicParticle*,
			 G4double tmin,
			 G4double maxEnergy) override;

  inline void SetLowestKineticEnergy(G4double e);

  G4double MinPrimaryEnergy(const G4Material*, const G4ParticleDefinition*, 
                            G4double) override;

  // hide assignment operator
  G4MuBremsstrahlungModel & 
    operator=(const  G4MuBremsstrahlungModel &right) = delete;
  G4MuBremsstrahlungModel(const  G4MuBremsstrahlungModel&) = delete;

protected:

  G4double ComputMuBremLoss(G4double Z, G4double tkin, G4double cut);
  
  G4double ComputeMicroscopicCrossSection(G4double tkin,
					  G4double Z,
					  G4double cut);

  virtual G4double ComputeDMicroscopicCrossSection(G4double tkin,
						   G4double Z,
						   G4double gammaEnergy);

  void SetParticle(const G4ParticleDefinition*);

protected:

  const G4ParticleDefinition* particle = nullptr;
  G4ParticleDefinition* theGamma = nullptr;
  G4ParticleChangeForLoss* fParticleChange = nullptr;
  G4NistManager* nist = nullptr;

  G4double mass = 1.0;
  G4double rmass = 1.0;
  G4double cc = 1.0;
  G4double coeff = 1.0;
  G4double sqrte;
  G4double bh = 202.4;
  G4double bh1 = 446.;
  G4double btf = 183.;
  G4double btf1 = 1429.;

  G4double lowestKinEnergy;
  G4double minThreshold;

  static const G4double xgi[6];
  static const G4double wgi[6];
  static G4double fDN[93];
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void G4MuBremsstrahlungModel::SetLowestKineticEnergy(G4double e) 
{
  lowestKinEnergy = e;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
