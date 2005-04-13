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
// $Id: G4KleinNishinaCompton.hh,v 1.3 2005-04-13 13:48:59 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4KleinNishinaCompton
//
// Author:        Vladimir Ivanchenko on base of Michel Maire code
//
// Creation date: 15.03.2005
//
// Modifications:
//
//
// Class Description:
//
// Implementation of gamma Compton scatteribg on free electron
// 

// -------------------------------------------------------------------
//

#ifndef G4KleinNishinaCompton_h
#define G4KleinNishinaCompton_h 1

#include "G4VEmModel.hh"

class G4ParticleChangeForLoss;

class G4KleinNishinaCompton : public G4VEmModel
{

public:

  G4KleinNishinaCompton(const G4ParticleDefinition* p = 0, 
			const G4String& nam = "Klein-Nishina");

  virtual ~G4KleinNishinaCompton();

  void Initialise(const G4ParticleDefinition*, const G4DataVector&);

  G4double ComputeCrossSectionPerAtom(
                                const G4ParticleDefinition*,
                                      G4double kinEnergy, 
                                      G4double Z, 
                                      G4double A, 
                                      G4double cut,
                                      G4double emax);

  std::vector<G4DynamicParticle*>* SampleSecondaries(
                                const G4MaterialCutsCouple*,
                                const G4DynamicParticle*,
                                      G4double tmin,
                                      G4double maxEnergy);

private:

  // hide assignment operator
  G4KleinNishinaCompton & operator=(const  G4KleinNishinaCompton &right);
  G4KleinNishinaCompton(const  G4KleinNishinaCompton&);

  G4ParticleDefinition*     theGamma;
  G4ParticleDefinition*     theElectron;
  G4ParticleChangeForLoss*  fParticleChange;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
