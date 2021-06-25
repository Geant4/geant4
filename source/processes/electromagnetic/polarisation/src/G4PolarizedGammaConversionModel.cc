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
// -------------------------------------------------------------------
//
// Geant4 Class file
//
// File name:     G4PolarizedGammaConversionModel
//
// Author:        Karim Laihem
//
// Class Description:
//   Implementation of gamma conversion to e+e- in the field of a nucleus
//   including polarization transfer

#include "G4PolarizedGammaConversionModel.hh"

#include "G4DataVector.hh"
#include "G4DynamicParticle.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4PolarizedGammaConversionXS.hh"
#include "G4PolarizationHelper.hh"
#include "G4StokesVector.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4PolarizedGammaConversionModel::G4PolarizedGammaConversionModel(
  const G4ParticleDefinition* pd, const G4String& nam)
  : G4BetheHeitlerModel(pd, nam)
  , fCrossSectionCalculator(nullptr)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4PolarizedGammaConversionModel::~G4PolarizedGammaConversionModel()
{
  if(fCrossSectionCalculator)
    delete fCrossSectionCalculator;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void G4PolarizedGammaConversionModel::Initialise(const G4ParticleDefinition* pd,
                                                 const G4DataVector& dv)
{
  G4BetheHeitlerModel::Initialise(pd, dv);
  if(!fCrossSectionCalculator)
    fCrossSectionCalculator = new G4PolarizedGammaConversionXS();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4PolarizedGammaConversionModel::SampleSecondaries(
  std::vector<G4DynamicParticle*>* vdp, const G4MaterialCutsCouple* couple,
  const G4DynamicParticle* dp, G4double tmin, G4double maxEnergy)
{
  G4BetheHeitlerModel::SampleSecondaries(vdp, couple, dp, tmin, maxEnergy);

  if(vdp && vdp->size() > 0)
  {
    G4double gamEnergy0 = dp->GetKineticEnergy();
    G4double lepEnergy1 = (*vdp)[0]->GetKineticEnergy();
    G4double sintheta =
      dp->GetMomentumDirection().cross((*vdp)[0]->GetMomentumDirection()).mag();
    if(sintheta > 1.)
      sintheta = 1.;

    G4StokesVector beamPol = G4StokesVector(dp->GetPolarization());
    beamPol.SetPhoton();

    // determine interaction plane
    G4ThreeVector nInteractionFrame = G4PolarizationHelper::GetFrame(
      dp->GetMomentumDirection(), (*vdp)[0]->GetMomentumDirection());

    // transform polarization into interaction frame
    beamPol.InvRotateAz(nInteractionFrame, dp->GetMomentumDirection());

    // calculate polarization transfer
    fCrossSectionCalculator->SetMaterial(
      GetCurrentElement()->GetN(),  // number of nucleons
      GetCurrentElement()->GetZ(), GetCurrentElement()->GetfCoulomb());
    fCrossSectionCalculator->Initialize(gamEnergy0, lepEnergy1, sintheta,
                                        beamPol, G4StokesVector::ZERO);

    // determine final state polarization
    G4StokesVector lep1Pol = fCrossSectionCalculator->GetPol2();
    lep1Pol.RotateAz(nInteractionFrame, (*vdp)[0]->GetMomentumDirection());
    (*vdp)[0]->SetPolarization(lep1Pol.p1(), lep1Pol.p2(), lep1Pol.p3());

    size_t num = vdp->size();
    if(num != 2)
    {
      G4ExceptionDescription ed;
      ed << " WARNING " << num
         << " secondaries in polarized pairproduction not supported!\n";
      G4Exception("G4PolarizedGammaConversionModel::SampleSecondaries",
                  "pol018", JustWarning, ed);
    }
    for(size_t i = 1; i < num; ++i)
    {
      G4StokesVector lep2Pol = fCrossSectionCalculator->GetPol3();
      lep2Pol.RotateAz(nInteractionFrame, (*vdp)[i]->GetMomentumDirection());
      (*vdp)[i]->SetPolarization(lep2Pol.p1(), lep2Pol.p2(), lep2Pol.p3());
    }
  }
}
