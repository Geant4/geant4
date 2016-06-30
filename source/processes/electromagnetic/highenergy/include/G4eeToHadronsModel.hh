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
// $Id: G4eeToHadronsModel.hh 97391 2016-06-02 10:08:45Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4eeToHadronsModel
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 25.10.2003
//
// Modifications:
// 08-04-05 Major optimisation of internal interfaces (V.Ivantchenko)
// 18-05-05 Use optimized interfaces (V.Ivantchenko)
//

//
// Class Description:
//

// -------------------------------------------------------------------
//

#ifndef G4eeToHadronsModel_h
#define G4eeToHadronsModel_h 1

#include "G4VEmModel.hh"

class G4PhysicsVector;
class G4Vee2hadrons;

class G4eeToHadronsModel : public G4VEmModel
{

public:

  explicit G4eeToHadronsModel(G4Vee2hadrons*, G4int ver=0,
                     const G4String& nam = "eeToHadrons");

  virtual ~G4eeToHadronsModel();

  virtual void Initialise(const G4ParticleDefinition*, 
                          const G4DataVector&) override;

  virtual G4double CrossSectionPerVolume(const G4Material*,
					 const G4ParticleDefinition*,
					 G4double kineticEnergy,
					 G4double cutEnergy,
					 G4double maxEnergy) override;

  virtual G4double ComputeCrossSectionPerAtom(
                                       const G4ParticleDefinition*,
                                       G4double kineticEnergy,
                                       G4double Z, G4double A,
                                       G4double cutEnergy = 0.0,
                                       G4double maxEnergy = DBL_MAX) override;

  virtual G4double ComputeCrossSectionPerElectron(
                                       const G4ParticleDefinition*,
                                       G4double kineticEnergy,
                                       G4double cutEnergy = 0.0,
                                       G4double maxEnergy = DBL_MAX);

  virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
				 const G4MaterialCutsCouple*,
				 const G4DynamicParticle*,
				 G4double tmin = 0.0,
				 G4double maxEnergy = DBL_MAX) override;

  G4DynamicParticle* GenerateCMPhoton(G4double);

  inline G4double PeakEnergy() const;

private:

  void ComputeCMCrossSectionPerElectron();

  // hide assignment operator
  G4eeToHadronsModel & operator=(const  G4eeToHadronsModel &right) = delete;
  G4eeToHadronsModel(const  G4eeToHadronsModel&) = delete;

  G4Vee2hadrons*        model;
  G4ParticleDefinition* theGamma;
  G4PhysicsVector*      crossPerElectron;
  G4PhysicsVector*      crossBornPerElectron;
  G4bool                isInitialised;
  G4int                 nbins;
  G4int                 verbose;

  G4double              lowKinEnergy;
  G4double              peakKinEnergy;
  G4double              highKinEnergy;

  G4double              emin;
  G4double              epeak;
  G4double              emax;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4eeToHadronsModel::PeakEnergy() const
{
  return peakKinEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
