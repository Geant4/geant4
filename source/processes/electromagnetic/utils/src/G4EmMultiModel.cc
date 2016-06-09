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
// $Id: G4EmMultiModel.cc,v 1.6 2007/05/22 17:31:58 vnivanch Exp $
// GEANT4 tag $Name: geant4-09-01 $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:   G4EmMultiModel
//
// Author:        Vladimir Ivanchenko
// 
// Creation date: 03.05.2004
//
// Modifications: 
// 15-04-05 optimize internal interface (V.Ivanchenko)
//

// Class Description:
//
// Energy loss model using several G4VEmModels

// -------------------------------------------------------------------
//


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4EmMultiModel.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4EmMultiModel::G4EmMultiModel(const G4String& nam)
  : G4VEmModel(nam),
  nModels(0)
{
  model.clear();
  tsecmin.clear();
  cross_section.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4EmMultiModel::~G4EmMultiModel()
{
  if(nModels) {
    for(G4int i=0; i<nModels; i++) {
      delete model[i];
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmMultiModel::Initialise(const G4ParticleDefinition* p, 
                                const G4DataVector& cuts)
{
  if(nModels) {
    for(G4int i=0; i<nModels; i++) {
      (model[i])->Initialise(p, cuts);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmMultiModel::MinEnergyCut(const G4ParticleDefinition* p,
                                      const G4MaterialCutsCouple* couple)
{
  G4double cut = DBL_MAX;
  if(nModels) {
    cut = (model[0])->MinEnergyCut(p, couple);
  } 
  return cut;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmMultiModel::ComputeDEDX(const G4MaterialCutsCouple* couple,
                                     const G4ParticleDefinition* p,
                                           G4double kineticEnergy,
                                           G4double cutEnergy)
{
  G4double dedx  = 0.0;

  if(nModels) {
    dedx =  (model[0])->ComputeDEDX(couple, p, cutEnergy, kineticEnergy);
  } 

  return dedx;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmMultiModel::CrossSection(const G4MaterialCutsCouple* couple,
                                      const G4ParticleDefinition* p,
                                            G4double kineticEnergy,
                                            G4double cutEnergy,
                                            G4double maxKinEnergy)
{
  G4double cross     = 0.0;
  G4double t1        = cutEnergy;
  G4double t2        = cutEnergy;
  if(nModels) {
    for(G4int i=0; i<nModels; i++) {
      t1 = std::max(t2, tsecmin[i]);
      t2 = std::min(maxKinEnergy, tsecmin[i+1]);
      cross += (model[i])->CrossSection(couple, p, kineticEnergy, t1, t2);
    }
  } 
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmMultiModel::SampleSecondaries(std::vector<G4DynamicParticle*>* vdp,
				       const G4MaterialCutsCouple* couple,
				       const G4DynamicParticle* dp,
				       G4double tmin,
				       G4double maxEnergy)
{
  if(nModels > 0) {
    G4int i;
    G4double cross = 0.0;
    G4double t1    = tmin;
    G4double t2    = tmin;
    for(i=0; i<nModels; i++) {
      t1 = std::max(t2, tsecmin[i]);
      t2 = std::min(maxEnergy, tsecmin[i+1]);
      cross += (model[i])->CrossSection(couple, dp->GetDefinition(), 
                                        dp->GetKineticEnergy(), t1, t2);
      cross_section[i] = cross;
    }

    cross *= G4UniformRand();
    t2 = tmin;

    for(i=0; i<nModels; i++) {
      t1 = std::max(t2, tsecmin[i]);
      t2 = std::min(maxEnergy, tsecmin[i+1]);
      if(cross <= cross_section[i]) {
        (model[i])->SampleSecondaries(vdp, couple, dp, t1, t2);
        break;
      }
    }
  } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmMultiModel:: MaxSecondaryEnergy(const G4ParticleDefinition*,
					           G4double kinEnergy)
{
  G4cout << "Warning! G4EmMultiModel::"
         << "MaxSecondaryEnergy(const G4ParticleDefinition*,G4double kinEnergy)"
         << " should not be used!" << G4endl;
  return kinEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmMultiModel::DefineForRegion(const G4Region* r)
{
  if(nModels) {
    for(G4int i=0; i<nModels; i++) {(model[i])->DefineForRegion(r);}
  } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmMultiModel::AddModel(G4VEmModel* p, G4double tmin, G4double tmax)
{
  if(tmin < tmax && 0.0 < tmin) {

    if(nModels == 0) {
      tsecmin.push_back(tmin);
      tsecmin.push_back(tmax);
      cross_section.push_back(0.0);
      model.push_back(p);
      nModels++;

    } else {
      G4int i, j;
      G4bool increment = false;
      for(i=0; i<nModels; i++) {

        if(tmin < tsecmin[i]) {
          G4double t2 = std::min(tsecmin[i+1],tmax);
          if(tmin < t2) {
            tsecmin.push_back(0.0);
            cross_section.push_back(0.0);
            model.push_back(0);
            for(j=nModels; j>i; j--) {
              model[j] = model[j-1];
              tsecmin[j+1] = tsecmin[j];
	    } 
            model[i] = p;
            tsecmin[i+1] = t2;
            tsecmin[i]   = tmin;
            increment = true;
	  }
   	} else if(i == nModels-1) {
          G4double t1 = std::min(tsecmin[i+1],tmin);
          G4double t2 = std::max(tsecmin[i+1],tmax);
          if(t1 < t2) {
            tsecmin.push_back(t2);
            cross_section.push_back(0.0);
            model.push_back(p);
            increment = true;
          }
	}
      }
      if(increment) nModels++;
    }
  } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

