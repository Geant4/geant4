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
// $Id: G4PolarizedPEEffectModel.cc 96114 2016-03-16 18:51:33Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4PolarizedPEEffectModel
//
// Author:        Andreas Schaelicke & Karim Laihem
//
// Creation date: 22.02.2007
//
// Modifications:
//
// Class Description:
//
// Implementation of Photo electric effect 
// including polarization transfer from circularly polarised gammas

// -------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef NOIONIZATIONAS

#include "G4PolarizedPEEffectModel.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"
#include "Randomize.hh"
#include "G4DataVector.hh"
#include "G4PhysicsLogVector.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4PolarizedPEEffectCrossSection.hh"
#include "G4PolarizationHelper.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


G4PolarizedPEEffectModel::G4PolarizedPEEffectModel(const G4ParticleDefinition*,
						   const G4String& nam)
  : G4PEEffectFluoModel(nam),
    crossSectionCalculator(nullptr),
    verboseLevel(0)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PolarizedPEEffectModel::~G4PolarizedPEEffectModel()
{
  if (crossSectionCalculator) delete crossSectionCalculator;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PolarizedPEEffectModel::Initialise(const G4ParticleDefinition* pd,
				     const G4DataVector& dv)
{
  G4PEEffectFluoModel::Initialise(pd,dv);
  if (!crossSectionCalculator)
    crossSectionCalculator = new G4PolarizedPEEffectCrossSection();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4PolarizedPEEffectModel::SampleSecondaries(std::vector<G4DynamicParticle*>* vdp,
						 const G4MaterialCutsCouple* couple,
						 const G4DynamicParticle* dp,
						 G4double tmin,
						 G4double maxEnergy)
{
  //  std::vector<G4DynamicParticle*>* vdp = 
  G4PEEffectFluoModel::SampleSecondaries(vdp,couple, dp, tmin, maxEnergy);

  if (verboseLevel >= 1) {
    G4cout << "G4PolarizedPEEffectModel::SampleSecondaries" << G4endl;
  }

  if(vdp && vdp->size()>0) {
    G4double gamEnergy0 = dp->GetKineticEnergy();
    G4double lepEnergy1 = (*vdp)[0]->GetKineticEnergy();
    G4double sintheta   = dp->GetMomentumDirection().cross((*vdp)[0]->GetMomentumDirection()).mag();
    if (sintheta>1.) sintheta=1.;

    G4StokesVector beamPol = dp->GetPolarization();
    beamPol.SetPhoton();
    //    G4cout<<" beamPol "<<beamPol<<G4endl;

    // determine interaction plane
    G4ThreeVector  nInteractionFrame = 
      G4PolarizationHelper::GetFrame(dp->GetMomentumDirection(), 
				     (*vdp)[0]->GetMomentumDirection());
    //    G4cout<<" nInteractionFrame = "<<nInteractionFrame<<G4endl;
    if (dp->GetMomentumDirection().cross((*vdp)[0]->GetMomentumDirection()).mag()<1.e-10) {
      //      G4cout<<" nInteractionFrame not well defined "<<G4endl;
      //      G4cout<<" choosing random interaction frame "<<G4endl;

      nInteractionFrame = G4PolarizationHelper::GetRandomFrame(dp->GetMomentumDirection());
      //      G4cout<<"new nInteractionFrame = "<<nInteractionFrame<<G4endl;
    }
    /*
    else {
      G4ThreeVector mom1=dp->GetMomentumDirection();
      G4ThreeVector mom2=(*vdp)[0]->GetMomentumDirection();
      G4cout<<" mom1 = "<<mom1<<G4endl;
      G4cout<<" mom2 = "<<mom2<<G4endl;
      G4ThreeVector x=mom1.cross(mom2);
      G4cout<<" mom1 x mom2 = "<<x<<" "<<x.mag()<<G4endl;
      G4cout<<"  norm       = "<<(1./x.mag()*x)<<" "<<G4endl;
    }
    */

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
    //    G4cout<<" lepPol "<<lep1Pol<<G4endl;
    lep1Pol.RotateAz(nInteractionFrame,(*vdp)[0]->GetMomentumDirection());
    (*vdp)[0]->SetPolarization(lep1Pol.p1(),
			       lep1Pol.p2(),
			       lep1Pol.p3());

    //    G4cout<<" lepPol "<<lep1Pol<<G4endl;
    size_t num = vdp->size();
    if (num!=1) G4cout<<" WARNING "<<num<<" secondaries in polarized photo electric effect not supported!\n"; 
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#endif //NOIONIZATIONAS
