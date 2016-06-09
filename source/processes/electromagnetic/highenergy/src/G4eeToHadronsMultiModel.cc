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
// $Id: G4eeToHadronsMultiModel.cc,v 1.2 2005/11/23 19:17:26 vnivanch Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4eeToHadronsMultiModel
//
// Author:        Vladimir Ivanchenko on base of Michel Maire code
//
// Creation date: 02.08.2004
//
// Modifications:
// 08-11-04 Migration to new interface of Store/Retrieve tables (V.Ivantchenko)
// 08-04-05 Major optimisation of internal interfaces (V.Ivantchenko)
//

//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4eeToHadronsMultiModel.hh"
#include "G4eeToTwoPiModel.hh"
#include "G4eeCrossSections.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

G4eeToHadronsMultiModel::G4eeToHadronsMultiModel(G4int ver, const G4String& name)
  : G4VEmModel(name),
    csFactor(1.0),
    nModels(0),
    verbose(ver),
    isInitialised(false)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4eeToHadronsMultiModel::~G4eeToHadronsMultiModel()
{
  G4int n = models.size();
  if(n>0) {
    for(G4int i=0; i<n; i++) {
      delete models[i];
    }
  }
  delete cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eeToHadronsMultiModel::Initialise(const G4ParticleDefinition* p, const G4DataVector& v)
{
  if(!isInitialised) {
    isInitialised = true;

    thKineticEnergy = DBL_MAX;
    maxKineticEnergy = HighEnergyLimit();

    cross = new G4eeCrossSections();
    G4eeToHadronsModel* model = new G4eeToHadronsModel(new G4eeToTwoPiModel(cross), verbose);
    models.push_back(model);
    model->SetHighEnergyLimit(maxKineticEnergy);
    model->Initialise(p, v);
    G4double emin = model->LowEnergyLimit();
    if(emin < thKineticEnergy) thKineticEnergy = emin;
    ekinMin.push_back(emin);
    ekinMax.push_back(model->HighEnergyLimit());
    ekinPeak.push_back(model->PeakEnergy());
    cumSum.push_back(0.0);
    nModels = 1;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eeToHadronsMultiModel::PrintInfo()
{
  if(verbose > 0) {
    G4cout << "      e+ annihilation into hadrons active above "
           << thKineticEnergy/GeV << " GeV"
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eeToHadronsMultiModel::SetCrossSecFactor(G4double fac)
{
  if(fac > 1.0) {
    csFactor = fac;
    if(verbose > 0)
      G4cout << "### G4eeToHadronsMultiModel: The cross section for G4eeToHadronsMultiModel is  "
             << "increased by the Factor= " << csFactor << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
