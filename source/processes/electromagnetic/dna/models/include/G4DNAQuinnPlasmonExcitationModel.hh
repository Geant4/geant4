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
// Created on 2016/04/08
//
// Authors: D. Sakata, S. Incerti
//
// This class perform transmission term of volume plasmon excitation,
// based on Quinn Model, see Phys. Rev. vol 126, number 4 (1962)

#ifndef G4DNAQuinnPlasmonExcitationModel_h
#define G4DNAQuinnPlasmonExcitationModel_h 1

#include "G4VEmModel.hh"
#include "G4Electron.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4ProductionCutsTable.hh"

class G4DNAQuinnPlasmonExcitationModel: public G4VEmModel
{

public:

  G4DNAQuinnPlasmonExcitationModel(const G4ParticleDefinition* p = 0,
                              const G4String& nam = "DNAQuinnPlasmonExcitationModel");

  virtual ~G4DNAQuinnPlasmonExcitationModel();

  virtual void Initialise(const G4ParticleDefinition*,
                          const G4DataVector& = *(new G4DataVector()));

  virtual G4double CrossSectionPerVolume(const G4Material* material,
                                         const G4ParticleDefinition* p,
                                         G4double ekin,
                                         G4double emin,
                                         G4double emax);

  virtual G4double GetCrossSection(const G4Material* material,
                                          const G4ParticleDefinition*,
                                          G4double kineticEnergy);

  virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
                                 const G4MaterialCutsCouple*,
                                 const G4DynamicParticle*,
                                 G4double tmin,
                                 G4double maxEnergy);
  
  inline void SelectStationary(G4bool input);

protected:

  G4ParticleChangeForGamma* fParticleChangeForGamma;

private:

  G4double fLowEnergyLimit=0.;
  G4double fHighEnergyLimit=0.;

  G4bool   isInitialised=false;
  G4bool   statCode=false;
  G4int    verboseLevel=0;
  G4int    nValenceElectron[100];

  const  std::vector<G4double>* fpMaterialDensity=nullptr;

  G4int GetNValenceElectron(G4int z);
   
  G4DNAQuinnPlasmonExcitationModel & operator=(const  G4DNAQuinnPlasmonExcitationModel &right);
  G4DNAQuinnPlasmonExcitationModel(const  G4DNAQuinnPlasmonExcitationModel&);

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4DNAQuinnPlasmonExcitationModel::SelectStationary(G4bool input)
{
  statCode = input;
}

#endif
