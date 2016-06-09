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
// $Id$
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4HeatedKleinNishinaCompton
//
// Author:        Vladimir Grichine on base of M. Maire and V. Ivanchenko code
//
// Creation date: 15.03.2009
//
// Modifications:
// 
//
//
// Class Description:
//
// Implementation of gamma Compton scattering on heated electrons 
// Maxwell distributed with the temperature fTemperature (default = 1*keV)
// 

// -------------------------------------------------------------------
//

#ifndef G4HeatedKleinNishinaCompton_h
#define G4HeatedKleinNishinaCompton_h 1

#include "G4VEmModel.hh"

class G4ParticleChangeForGamma;

class G4HeatedKleinNishinaCompton : public G4VEmModel
{

public:

  G4HeatedKleinNishinaCompton(const G4ParticleDefinition* p = 0, 
			const G4String& nam = "Heated-Klein-Nishina");

  virtual ~G4HeatedKleinNishinaCompton();

  virtual void Initialise(const G4ParticleDefinition*, const G4DataVector&);

  virtual G4double ComputeCrossSectionPerAtom(
                                const G4ParticleDefinition*,
                                      G4double kinEnergy, 
                                      G4double Z, 
                                      G4double A, 
                                      G4double cut,
                                      G4double emax);

  virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
				 const G4MaterialCutsCouple*,
				 const G4DynamicParticle*,
				 G4double tmin,
				 G4double maxEnergy);

  void     SetElectronTemperature(G4double t){fTemperature=t;};
  G4double GetElectronTemperature(){return fTemperature;};

protected:

  G4ParticleDefinition*     theGamma;
  G4ParticleDefinition*     theElectron;
  G4ParticleChangeForGamma* fParticleChange;
  G4double                  lowestGammaEnergy;
  G4double                  fTemperature;  // electron temperature in energy units

private:

  // hide assignment operator
  G4HeatedKleinNishinaCompton & operator=(const  G4HeatedKleinNishinaCompton &right);
  G4HeatedKleinNishinaCompton(const  G4HeatedKleinNishinaCompton&);

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
