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
// $Id: G4PolarizedGammaConversionModel.cc 96114 2016-03-16 18:51:33Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4PolarizedGammaConversionModel
//
// Author:        Karim Laihem
//
// Creation date: 19.04.2005
//
// Modifications:
// 21-08-06 Modified to fit in g4.8.1 framework (A.Schaelicke)
// 19-03-07 Add initialisation of crossSectionCalculator (VI)
//
// Class Description:
//
// Implementation of gamma convertion to e+e- in the field of a nucleus 
// including polarization transfer

// -------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4PolarizedGammaConversionModel.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"
#include "Randomize.hh"
#include "G4DataVector.hh"
#include "G4PhysicsLogVector.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4PolarizedPairProductionCrossSection.hh"
#include "G4PolarizationHelper.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


G4PolarizedGammaConversionModel::G4PolarizedGammaConversionModel(const G4ParticleDefinition* pd,
					 const G4String& nam)
  : G4BetheHeitlerModel(pd,nam), crossSectionCalculator(nullptr)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PolarizedGammaConversionModel::~G4PolarizedGammaConversionModel()
{
  if (crossSectionCalculator) delete crossSectionCalculator;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PolarizedGammaConversionModel::Initialise(const G4ParticleDefinition* pd,
				     const G4DataVector& dv)
{
  G4BetheHeitlerModel::Initialise(pd,dv);
  if (!crossSectionCalculator)
    crossSectionCalculator = new G4PolarizedPairProductionCrossSection();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4PolarizedGammaConversionModel::SampleSecondaries(std::vector<G4DynamicParticle*>* vdp,
							const G4MaterialCutsCouple* couple,
							const G4DynamicParticle* dp,
							G4double tmin,
							G4double maxEnergy)
{
  G4BetheHeitlerModel::SampleSecondaries(vdp, couple, dp, tmin, maxEnergy);
 
  if(vdp && vdp->size()>0) {
    G4double gamEnergy0 = dp->GetKineticEnergy();
    G4double lepEnergy1 = (*vdp)[0]->GetKineticEnergy();
    G4double sintheta   = dp->GetMomentumDirection().cross((*vdp)[0]->GetMomentumDirection()).mag();
    if (sintheta>1.) sintheta=1.;

    G4StokesVector beamPol = dp->GetPolarization();
    beamPol.SetPhoton();

    // determine interaction plane
    G4ThreeVector  nInteractionFrame = 
      G4PolarizationHelper::GetFrame(dp->GetMomentumDirection(), 
				     (*vdp)[0]->GetMomentumDirection());

    // transform polarization into interaction frame
    beamPol.InvRotateAz(nInteractionFrame,dp->GetMomentumDirection());

    // calulcate polarization transfer
    crossSectionCalculator->SetMaterial(GetCurrentElement()->GetN(), // number of nucleons
					GetCurrentElement()->GetZ(), 
					GetCurrentElement()->GetfCoulomb());
    crossSectionCalculator->Initialize(gamEnergy0, lepEnergy1, sintheta,
				       beamPol, G4StokesVector::ZERO);

    // deterimine final state polarization
    G4StokesVector lep1Pol = crossSectionCalculator->GetPol2();
    lep1Pol.RotateAz(nInteractionFrame,(*vdp)[0]->GetMomentumDirection());
    (*vdp)[0]->SetPolarization(lep1Pol.p1(),
			       lep1Pol.p2(),
			       lep1Pol.p3());

    size_t num = vdp->size();
    if (num!=2) G4cout<<" WARNING "<<num<<" secondaries in polarized pairproduction not supported!\n"; 
    for (size_t i =1; i<num; ++i) {
      G4StokesVector lep2Pol = crossSectionCalculator->GetPol3();
      lep2Pol.RotateAz(nInteractionFrame,(*vdp)[i]->GetMomentumDirection());
      (*vdp)[i]->SetPolarization(lep2Pol.p1(),
				 lep2Pol.p2(),
				 lep2Pol.p3());

    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
