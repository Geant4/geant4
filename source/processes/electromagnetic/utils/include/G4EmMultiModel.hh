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
// $Id: G4EmMultiModel.hh,v 1.1 2004/05/03 13:13:07 vnivanch Exp $
// GEANT4 tag $Name: geant4-06-02 $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4EmMultiModel
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 03.05.2004
//
// Modifications:
//
//
//
// Class Description:
//
// Energy loss model using several G4VEmModels

// -------------------------------------------------------------------
//

#ifndef G4EmMultiModel_h
#define G4EmMultiModel_h 1

#include "globals.hh"
#include "G4VEmModel.hh"
#include <vector>

class G4Region;
class G4PhysicsTable;
class G4DynamicParticle;

class G4EmMultiModel :  public G4VEmModel
{

public:

  G4EmMultiModel(const G4String& nam = "MultiModel");

  ~G4EmMultiModel();

  void Initialise(const G4ParticleDefinition*, const G4DataVector&);

  G4double HighEnergyLimit(const G4ParticleDefinition*) {return highKinEnergy;};

  G4double LowEnergyLimit(const G4ParticleDefinition*) {return lowKinEnergy;};

  void SetHighEnergyLimit(G4double e) {highKinEnergy = e;};

  void SetLowEnergyLimit(G4double e) {lowKinEnergy = e;};

  G4double MinEnergyCut(const G4ParticleDefinition*,
                        const G4MaterialCutsCouple*);

  G4bool IsInCharge(const G4ParticleDefinition*);

  G4double ComputeDEDX(const G4MaterialCutsCouple*,
                               const G4ParticleDefinition*,
                                     G4double kineticEnergy,
                                     G4double cutEnergy);

  G4double CrossSection(const G4MaterialCutsCouple*,
                                const G4ParticleDefinition*,
                                      G4double kineticEnergy,
                                      G4double cutEnergy,
                                      G4double maxEnergy);

  G4DynamicParticle* SampleSecondary(
                                const G4MaterialCutsCouple*,
                                const G4DynamicParticle*,
                                      G4double tmin,
                                      G4double tmax);

  std::vector<G4DynamicParticle*>* SampleSecondaries(
                                const G4MaterialCutsCouple*,
                                const G4DynamicParticle*,
                                      G4double tmin,
                                      G4double tmax);

  G4double MaxSecondaryEnergy(const G4DynamicParticle* dynParticle);

  void DefineForRegion(const G4Region*);

  void SetDynamicParticle(const G4DynamicParticle*);

  void AddModel(G4VEmModel*, G4double tmin, G4double tmax);

protected:

  G4double MaxSecondaryEnergy(const G4ParticleDefinition*,
    				    G4double kineticEnergy);

private: 

  //  hide assignment operator
  G4EmMultiModel & operator=(const  G4EmMultiModel &right);
  G4EmMultiModel(const  G4EmMultiModel&);

  G4int                         nModels;
  std::vector<G4VEmModel*>      model;
  G4DataVector                  tsecmin;
  G4DataVector                  cross_section;

  G4double                      highKinEnergy;
  G4double                      lowKinEnergy;

};

#endif

