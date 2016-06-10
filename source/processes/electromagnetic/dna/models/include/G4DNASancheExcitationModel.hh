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
// $Id: G4DNASancheExcitationModel.hh 93616 2015-10-27 08:59:17Z gcosmo $
// GEANT4 tag $Name:  $
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Created by Z. Francis

#ifndef G4DNASancheExcitationModel_h
#define G4DNASancheExcitationModel_h 1

#include <deque>
#include <CLHEP/Units/SystemOfUnits.h>

#include "G4VEmModel.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4Electron.hh"
#include "G4NistManager.hh"

class G4DNASancheExcitationModel : public G4VEmModel
{

public:

  G4DNASancheExcitationModel(const G4ParticleDefinition* p = 0,
                              const G4String& nam = "DNASancheExcitationModel");

  virtual ~G4DNASancheExcitationModel();

  virtual void Initialise(const G4ParticleDefinition*, const G4DataVector&);

  virtual G4double CrossSectionPerVolume(const G4Material* material,
                                         const G4ParticleDefinition* p,
                                         G4double ekin,
                                         G4double emin,
                                         G4double emax);

  virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
                                 const G4MaterialCutsCouple*,
                                 const G4DynamicParticle*,
                                 G4double tmin,
                                 G4double maxEnergy);

  // Cross section

  G4double PartialCrossSection(G4double energy, G4int level);
  G4double TotalCrossSection(G4double t);

  inline void ExtendLowEnergyLimit(G4double /*threshold*/);

  inline void SetVerboseLevel(int verbose)
  {
    verboseLevel = verbose;
  }

protected:

  G4ParticleChangeForGamma* fParticleChangeForGamma;

private:
  // Water density table
  const std::vector<G4double>* fpWaterDensity;

  G4double lowEnergyLimit;
  G4double highEnergyLimit;
  G4bool isInitialised;
  G4int verboseLevel;

  // Cross section

  G4int RandomSelect(G4double energy);
  G4int nLevels;
  G4double VibrationEnergy(G4int level);
  G4double Sum(G4double k);
  G4double LinInterpolate(G4double e1,
                          G4double e2,
                          G4double e,
                          G4double xs1,
                          G4double xs2);

  //
//  typedef std::map<double, std::map<double, double> > TriDimensionMap;
//  TriDimensionMap map1;
  std::vector<double> tdummyVec;
  std::vector<std::vector<double>> fEnergyLevelXS;
  std::vector<double> fEnergyTotalXS;

  //
  G4DNASancheExcitationModel & operator=(const G4DNASancheExcitationModel &right);
  G4DNASancheExcitationModel(const G4DNASancheExcitationModel&);

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void G4DNASancheExcitationModel::ExtendLowEnergyLimit(G4double threshold)
{
  lowEnergyLimit = threshold;
  if(lowEnergyLimit < 2 * CLHEP::eV)
    G4Exception("*** WARNING : the G4DNASancheExcitationModel class is not "
                "validated below 2 eV !",
                "", JustWarning, "");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
