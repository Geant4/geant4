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
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4LindhardSorensenIonModel
//
// Author:        Alexander Bagulya & Vladimir Ivanchenko
//
// Creation date: 16.04.2018
//
//
// Class Description:
//
// Implementation of ion ionisation energy loss and delta-electron 
// production by heavy charged particles according to
// J. Lindhard & A.H. Sorensen, Phys. Rev. A 53 (1996) 2443-2455

// -------------------------------------------------------------------
//

#ifndef G4LindhardSorensenIonModel_h
#define G4LindhardSorensenIonModel_h 1

#include <vector>

#include "G4VEmModel.hh"
#include "G4NistManager.hh"

class G4EmCorrections;
class G4ParticleChangeForLoss;
class G4LindhardSorensenData;
class G4BraggIonModel;
class G4BetheBlochModel;
class G4IonICRU73Data;

class G4LindhardSorensenIonModel : public G4VEmModel
{
public:

  explicit G4LindhardSorensenIonModel(const G4ParticleDefinition* p = nullptr,
				      const G4String& nam = "LindhardSorensen");

  ~G4LindhardSorensenIonModel() override;

  void Initialise(const G4ParticleDefinition*, const G4DataVector&) override;

  G4double MinEnergyCut(const G4ParticleDefinition*,
                        const G4MaterialCutsCouple* couple) override;

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
                                G4double cutEnergy) override;

  G4double GetChargeSquareRatio(const G4ParticleDefinition* p,
                                const G4Material* mat,
                                G4double kineticEnergy) override;

  G4double GetParticleCharge(const G4ParticleDefinition* p,
			     const G4Material* mat,
                             G4double kineticEnergy) override;

  void CorrectionsAlongStep(const G4MaterialCutsCouple* couple,
                            const G4DynamicParticle* dp,
                            const G4double& length,
                            G4double& eloss) override;

  void SampleSecondaries(std::vector<G4DynamicParticle*>*,
                         const G4MaterialCutsCouple*,
		         const G4DynamicParticle*,
			 G4double tmin,
	                 G4double maxEnergy) override;

  // hide assignment operator
  G4LindhardSorensenIonModel & operator=
  (const  G4LindhardSorensenIonModel &right) = delete;
  G4LindhardSorensenIonModel(const  G4LindhardSorensenIonModel&) = delete;

protected:

  G4double MaxSecondaryEnergy(const G4ParticleDefinition*,
                              G4double kinEnergy) override;

  inline G4double GetChargeSquareRatio() const;

  inline void SetChargeSquareRatio(G4double val);

private:

  void SetupParameters();

  void InitialiseLS();

  G4double ComputeDEDXPerVolumeLS(const G4Material*,
                                  const G4ParticleDefinition*,
				  G4double kinEnergy, G4double cutEnergy);

  inline void SetParticle(const G4ParticleDefinition* p);

  static const G4int MAXZION = 93;

  static G4IonICRU73Data* fIonData;
  static G4LindhardSorensenData* lsdata;
  static std::vector<G4float>*   fact[MAXZION];

  const G4ParticleDefinition* particle;
  G4ParticleDefinition*       theElectron;
  G4EmCorrections*            corr;
  G4ParticleChangeForLoss*    fParticleChange;
  G4NistManager*              nist;
  G4BraggIonModel*            fBraggModel;
  G4BetheBlochModel*          fBBModel;

  G4int Zin = 1;
  G4double mass = 0.0;
  G4double tlimit = DBL_MAX;
  G4double spin = 0.0;
  G4double magMoment2 = 0.0;
  G4double chargeSquare = 1.0;
  G4double charge = 1.0;
  G4double eRatio = 0.0;
  G4double pRatio = 1.0;
  G4double formfact = 0.0;
  G4double twoln10;
  G4double fElimit;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void 
G4LindhardSorensenIonModel::SetParticle(const G4ParticleDefinition* p)
{
  if(particle != p) {
    particle = p;
    SetupParameters();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4LindhardSorensenIonModel::GetChargeSquareRatio() const
{
  return chargeSquare;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4LindhardSorensenIonModel::SetChargeSquareRatio(G4double val)
{
  chargeSquare = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
