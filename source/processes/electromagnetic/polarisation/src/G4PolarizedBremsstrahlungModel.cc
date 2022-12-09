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
// Geant4 class file
//
// File name:     G4PolarizedBremsstrahlungModel
//
// Author:        Karim Laihem

#include "G4PolarizedBremsstrahlungModel.hh"

#include "G4ParticleChangeForLoss.hh"
#include "G4PolarizationHelper.hh"
#include "G4PolarizedBremsstrahlungXS.hh"

G4PolarizedBremsstrahlungModel::G4PolarizedBremsstrahlungModel(
  const G4ParticleDefinition* p, const G4String& nam)
  : G4SeltzerBergerModel(p, nam)
  , fCrossSectionCalculator(nullptr)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4PolarizedBremsstrahlungModel::~G4PolarizedBremsstrahlungModel()
{
  if(fCrossSectionCalculator)
    delete fCrossSectionCalculator;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void G4PolarizedBremsstrahlungModel::Initialise(const G4ParticleDefinition* p,
                                                const G4DataVector& d)
{
  G4SeltzerBergerModel::Initialise(p, d);
  if(!fCrossSectionCalculator)
    fCrossSectionCalculator = new G4PolarizedBremsstrahlungXS();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// The emitted gamma energy is sampled using a parametrized formula
void G4PolarizedBremsstrahlungModel::SampleSecondaries(
  std::vector<G4DynamicParticle*>* vdp, const G4MaterialCutsCouple* couple,
  const G4DynamicParticle* dp, G4double tmin, G4double maxEnergy)
{
  G4SeltzerBergerModel::SampleSecondaries(vdp, couple, dp, tmin, maxEnergy);
  std::size_t num = vdp->size();

  if(num > 0)
  {
    G4double lepEnergy0 = dp->GetKineticEnergy();
    G4double gamEnergy1 = (*vdp)[0]->GetKineticEnergy();
    G4double sintheta =
      dp->GetMomentumDirection().cross((*vdp)[0]->GetMomentumDirection()).mag();
    if(sintheta > 1.)
      sintheta = 1.;

    G4StokesVector beamPol = G4StokesVector(dp->GetPolarization());

    // determine interaction plane
    G4ThreeVector nInteractionFrame = G4PolarizationHelper::GetFrame(
      dp->GetMomentumDirection(),
      fParticleChange->GetProposedMomentumDirection());

    // transform polarization into interaction frame
    beamPol.InvRotateAz(nInteractionFrame, dp->GetMomentumDirection());

    // calculate polarization transfer
    fCrossSectionCalculator->SetMaterial(
      GetCurrentElement()->GetN(),  // number of nucleons
      GetCurrentElement()->GetZ(), GetCurrentElement()->GetfCoulomb());
    fCrossSectionCalculator->Initialize(lepEnergy0, gamEnergy1, sintheta,
                                        beamPol, G4StokesVector::ZERO);

    // determine final state polarization
    G4StokesVector newBeamPol = fCrossSectionCalculator->GetPol2();
    newBeamPol.RotateAz(nInteractionFrame,
                        fParticleChange->GetProposedMomentumDirection());
    fParticleChange->ProposePolarization(newBeamPol);

    if(num != 1)
    {
      G4ExceptionDescription ed;
      ed << num << " secondaries in polarized bremsstrahlung not supported!\n";
      G4Exception("G4PolarizedBremsstrahlungModel::SampleSecondaries", "pol001",
                  JustWarning, ed);
    }
    for(std::size_t i = 0; i < num; ++i)
    {
      G4StokesVector photonPol = fCrossSectionCalculator->GetPol3();
      photonPol.SetPhoton();
      photonPol.RotateAz(nInteractionFrame, (*vdp)[i]->GetMomentumDirection());
      (*vdp)[i]->SetPolarization(photonPol.p1(), photonPol.p2(),
                                 photonPol.p3());
    }
  }
  return;
}
