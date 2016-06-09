//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: G4eeToHadronsModel.hh,v 1.3 2005/05/18 10:12:32 vnivanch Exp $
// GEANT4 tag $Name: geant4-07-01 $
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

  G4eeToHadronsModel(const G4Vee2hadrons*, G4int ver=0,
                     const G4String& nam = "eeToHadrons");

  virtual ~G4eeToHadronsModel();

  virtual void Initialise(const G4ParticleDefinition*, const G4DataVector&);

  G4double PeakEnergy() const;

  virtual G4double CrossSectionPerVolume(const G4Material*,
                                         const G4ParticleDefinition*,
                                         G4double kineticEnergy,
                                         G4double cutEnergy = 0.0,
                                         G4double maxEnergy = DBL_MAX);

  virtual std::vector<G4DynamicParticle*>* SampleSecondaries(
                                const G4MaterialCutsCouple*,
                                const G4DynamicParticle*,
                                      G4double tmin = 0.0,
                                      G4double maxEnergy = DBL_MAX);

  G4DynamicParticle* GenerateCMPhoton(G4double);

private:

  void ComputeCMCrossSectionPerElectron();

  // hide assignment operator
  G4eeToHadronsModel & operator=(const  G4eeToHadronsModel &right);
  G4eeToHadronsModel(const  G4eeToHadronsModel&);

  const G4Vee2hadrons*  model;
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
