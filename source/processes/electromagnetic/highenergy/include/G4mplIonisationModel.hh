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
// $Id: G4mplIonisationModel.hh,v 1.1 2006-10-25 17:37:44 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4mplIonisationModel
//
// Author:        Vladimir Ivanchenko 
//
// Creation date: 06.09.2005
//
// Modifications:
//
//
// Class Description:
//
// Implementation of model of energy loss of the magnetic monopole

// -------------------------------------------------------------------
//

#ifndef G4mplIonisationModel_h
#define G4mplIonisationModel_h 1

#include "G4VEmModel.hh"

class G4Monopole;
class G4ParticleChangeForLoss;

class G4mplIonisationModel : public G4VEmModel
{

public:

  G4mplIonisationModel(const G4ParticleDefinition* p = 0, const G4String& nam = "mplIonisation");

  virtual ~G4mplIonisationModel();

  virtual void Initialise(const G4ParticleDefinition*, const G4DataVector&);

  virtual G4double ComputeDEDXPerVolume(const G4Material*,
					const G4ParticleDefinition*,
					G4double kineticEnergy,
					G4double cutEnergy);

  virtual std::vector<G4DynamicParticle*>* SampleSecondaries(
                                const G4MaterialCutsCouple*,
                                const G4DynamicParticle*,
                                      G4double tmin,
                                      G4double maxEnergy);


private:

  void SetParticle(const G4ParticleDefinition* p);

  // hide assignment operator
  G4mplIonisationModel & operator=(const  G4mplIonisationModel &right);
  G4mplIonisationModel(const  G4mplIonisationModel&);

  const G4Monopole*           monopole;
  G4ParticleChangeForLoss*    fParticleChange;

  G4double mass;
  G4double magCharge;
  G4double chargeSquare;
  G4double twoln10;
  G4double beta2low;
  G4double beta2lim;
  G4double bg2lim;
  G4double factlow;
  G4int    nmpl;
};

#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline std::vector<G4DynamicParticle*>* G4mplIonisationModel::SampleSecondaries(
                             const G4MaterialCutsCouple*,
                             const G4DynamicParticle*,
                                   G4double,
                                   G4double)
{
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
