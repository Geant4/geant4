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
// $Id: G4PolarizedMollerBhabhaModel.cc 96114 2016-03-16 18:51:33Z gcosmo $
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
// File name:     G4PolarizedMollerBhabhaModel
//
// Author:        A.Schaelicke on base of Vladimir Ivanchenko code
//
// Creation date: 10.11.2005
//
// Modifications:
//
// 20-08-05, modified interface (A.Schaelicke)
//
// Class Description:
//
// Implementation of energy loss and delta-electron production by e+/e-
// (including polarization effects)
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4PolarizedMollerBhabhaModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4ParticleChangeForLoss.hh"
#include "Randomize.hh"

#include "G4PolarizationManager.hh"
#include "G4PolarizationHelper.hh"
#include "G4PolarizedBhabhaCrossSection.hh"
#include "G4PolarizedMollerCrossSection.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PolarizedMollerBhabhaModel::G4PolarizedMollerBhabhaModel(const G4ParticleDefinition* p,
                                         const G4String& nam)
  : G4MollerBhabhaModel(p,nam)
{

  //   G4cout<<" particle==electron "<<(p==theElectron)<<G4endl;
  isElectron=(p==theElectron);  // necessary due to wrong order in G4MollerBhabhaModel constructor!

  if (p==nullptr) { 
    
  }
  if (!isElectron) {
    G4cout<<" buildBhabha cross section "<<isElectron<<G4endl;
    crossSectionCalculator = new G4PolarizedBhabhaCrossSection();
  }  else {
    G4cout<<" buildMoller cross section "<<isElectron<<G4endl;
    crossSectionCalculator = new G4PolarizedMollerCrossSection();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4PolarizedMollerBhabhaModel::~G4PolarizedMollerBhabhaModel()
{
  if (crossSectionCalculator) {
    delete crossSectionCalculator;
  }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4PolarizedMollerBhabhaModel::ComputeCrossSectionPerElectron(
                                const G4ParticleDefinition* pd,
                                      G4double kinEnergy, 
                                      G4double cut,
                                      G4double emax)
{
  G4double xs = 
    G4MollerBhabhaModel::ComputeCrossSectionPerElectron(pd,kinEnergy,
							cut,emax);
//   G4cout<<"calc eIoni xsec "<<xs<<G4endl;
//   G4cout<<" "<<kinEnergy<<" "<<cut<<" "<<emax<<G4endl;
  G4double factor=1.;
  if (xs!=0.) {
    //    G4cout<<"calc asym"<<G4endl;
    G4double tmax = MaxSecondaryEnergy(pd, kinEnergy);
    tmax = std::min(emax, tmax);

    if (std::fabs(cut/emax-1.)<1.e-10) return xs;

    if(cut < tmax) {

      G4double xmin  = cut/kinEnergy;
      G4double xmax  = tmax/kinEnergy;
//       G4cout<<"calc asym "<<xmin<<","<<xmax<<G4endl;
      G4double gam   = kinEnergy/electron_mass_c2 + 1.0;

      G4double crossPol=crossSectionCalculator->
	TotalXSection(xmin,xmax,gam,
		      theBeamPolarization,
		      theTargetPolarization);
      G4double crossUnpol=crossSectionCalculator->
	TotalXSection(xmin,xmax,gam,
		      G4StokesVector::ZERO,
		      G4StokesVector::ZERO);
      if (crossUnpol>0.)  factor=crossPol/crossUnpol;
      //     G4cout<<" factor="<<factor<<G4endl;
    }
  }
  return xs*factor;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PolarizedMollerBhabhaModel::SampleSecondaries(std::vector<G4DynamicParticle*>* vdp,
						     const G4MaterialCutsCouple* ,
						     const G4DynamicParticle* dp,
						     G4double tmin,
						     G4double maxEnergy)
{
  // *** obtain and save target and beam polarization ***
  G4PolarizationManager * polarizationManger = G4PolarizationManager::GetInstance();

  const G4Track * aTrack = fParticleChange->GetCurrentTrack();

  // obtain polarization of the beam
  theBeamPolarization = dp->GetPolarization();

  // obtain polarization of the media
  G4VPhysicalVolume*  aPVolume  = aTrack->GetVolume();
  G4LogicalVolume*    aLVolume  = aPVolume->GetLogicalVolume();
  const G4bool targetIsPolarized = polarizationManger->IsPolarized(aLVolume);
  theTargetPolarization = polarizationManger->GetVolumePolarization(aLVolume);

  // transfer target polarization in interaction frame
  if (targetIsPolarized)
      theTargetPolarization.rotateUz(dp->GetMomentumDirection());




  G4double tmax = std::min(maxEnergy, MaxSecondaryKinEnergy(dp));
  if(tmin >= tmax) return;
  //  if(tmin > tmax) tmin = tmax;

  G4double polL = theBeamPolarization.z()*theTargetPolarization.z();
  		polL=std::fabs(polL);
  G4double polT = theBeamPolarization.x()*theTargetPolarization.x() +
  		  theBeamPolarization.y()*theTargetPolarization.y();
  		polT=std::fabs(polT);

  G4double kineticEnergy = dp->GetKineticEnergy();
  G4double energy = kineticEnergy + electron_mass_c2;
  G4double totalMomentum = std::sqrt(kineticEnergy*(energy + electron_mass_c2));
  G4double xmin   = tmin/kineticEnergy;
  G4double xmax   = tmax/kineticEnergy;
  G4double gam    = energy/electron_mass_c2;
  G4double gamma2 = gam*gam;
    G4double gmo  = gam - 1.;
    G4double gmo2 = gmo*gmo;
    G4double gmo3 = gmo2*gmo;
    G4double gpo  = gam + 1.;
    G4double gpo2 = gpo*gpo;
    G4double gpo3 = gpo2*gpo;
  G4double x, y, q, grej, grej2;
  G4double z = 0.;
  G4double xs = 0., phi =0.;
  G4ThreeVector direction = dp->GetMomentumDirection();

  //(Polarized) Moller (e-e-) scattering
  if (isElectron) {
    // *** dice according to polarized cross section
    G4double G = ((2.0*gam - 1.0)/gamma2)*(1. - polT - polL*gam);
    G4double H =  (sqr(gam - 1.0)/gamma2)*(1. + polT + polL*((gam + 3.)/(gam - 1.)));

    y  = 1.0 - xmax;
    grej  = 1.0 - G*xmax + xmax*xmax*(H + (1.0 - G*y)/(y*y));
    grej2 = 1.0 - G*xmin + xmin*xmin*(H + (1.0 - G*y)/(y*y));
    if (grej2 > grej) grej = grej2;
    G4double prefM = gamma2*classic_electr_radius*classic_electr_radius/(gmo2*(gam + 1.0));
    grej *= prefM;
    do {
      q = G4UniformRand();
      x = xmin*xmax/(xmin*(1.0 - q) + xmax*q);
      if (crossSectionCalculator) {
	crossSectionCalculator->Initialize(x,gam,phi,theBeamPolarization,
					   theTargetPolarization,1);
	xs=crossSectionCalculator->XSection(G4StokesVector::ZERO,
						     G4StokesVector::ZERO);
	z=xs*sqr(x)*4.;
	if (grej < z) {
	  G4cout<<"WARNING : error in Moller rejection routine! \n"
		<<" z = "<<z<<" grej="<<grej<<"\n";
	}
      } else {
	G4cout<<"No calculator in Moller scattering"<<G4endl;
      }
      // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
    } while(grej * G4UniformRand() > z);
    //Bhabha (e+e-) scattering
  } else {
    // *** dice according to polarized cross section
    y     = xmax*xmax;
    grej = 0.;
    grej += y*y*gmo3*(1. + (polL + polT)*(gam + 3.)/gmo);
    grej += -2.*xmin*xmin*xmin*gam*gmo2*(1. - (polL + polT)*(gam + 3.)/gmo);
    grej += y*y*gmo*(3.*gamma2 + 6.*gam + 4.)*(1. + (polL*(3.*gam + 1.)*(gamma2 + gam + 1.) + polT*((gam + 2.)*gamma2 + 1.))/(gmo*(3.*gam*(gam + 2.) + 4.)));
    grej /= gpo3;
    grej += -xmin*(2.*gamma2 + 4.*gam + 1.)*(1. - gam*(polL*(2.*gam + 1.) + polT)/(2.*gam*(gam + 2.) + 1.))/gpo2;
    grej += gamma2/(gamma2 - 1.);
    G4double prefB = classic_electr_radius*classic_electr_radius/(gam - 1.0);
    grej *= prefB;

    do {
      q  = G4UniformRand();
      x  = xmin*xmax/(xmin*(1.0 - q) + xmax*q);
      if (crossSectionCalculator) {
	crossSectionCalculator->Initialize(x,gam,phi,theBeamPolarization,
					   theTargetPolarization,1);
	xs=crossSectionCalculator->XSection(G4StokesVector::ZERO,
						     G4StokesVector::ZERO);
	z=xs*sqr(x)*4.;
      } else {
	G4cout<<"No calculator in Bhabha scattering"<<G4endl;
      }

      if(z > grej) {
	G4cout<<"&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"<<G4endl;
        G4cout << "G4PolarizedMollerBhabhaModel::SampleSecondaries Warning! "<<G4endl
               << "Majorant " << grej << " < "
               << z << " for x= " << x<<G4endl
               << " e+e- (Bhabha) scattering"<<" at KinEnergy "<<kineticEnergy<<G4endl;
	G4cout<<"&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"<<G4endl;
      }
      // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
    } while(grej * G4UniformRand() > z);
  }
  //
  //
  // polar asymmetries (due to transverse polarizations)
  //
  //
  if (crossSectionCalculator) {
   // grej*=1./(sqr(x)*sqr(gamma2-1))*sqr(gam*(1+gam));
    grej=xs*2.;
    do {
      phi = twopi * G4UniformRand() ;
      crossSectionCalculator->Initialize(x,gam,phi,theBeamPolarization,
					           theTargetPolarization,1);
      xs=crossSectionCalculator->XSection(G4StokesVector::ZERO,
					  G4StokesVector::ZERO);
      if(xs > grej) {
	if (isElectron){
	  G4cout << "G4PolarizedMollerBhabhaModel::SampleSecondaries Warning! "<<G4endl
		 << "Majorant " << grej << " < "
		 << xs << " for phi= " << phi<<G4endl
		 << " e-e- (Moller) scattering"<< G4endl
		 <<"PHI DICING"<<G4endl;
	} else {
	  G4cout << "G4PolarizedMollerBhabhaModel::SampleSecondaries Warning! "<<G4endl
		 << "Majorant " << grej << " < "
		 << xs << " for phi= " << phi<<G4endl
		 << " e+e- (Bhabha) scattering"<< G4endl
		 <<"PHI DICING"<<G4endl;
	}
      }
      // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
    } while(grej * G4UniformRand() > xs);
  }

  // fix kinematics of delta electron
  G4double deltaKinEnergy = x * kineticEnergy;
  G4double deltaMomentum =
           std::sqrt(deltaKinEnergy * (deltaKinEnergy + 2.0*electron_mass_c2));
  G4double cost = deltaKinEnergy * (energy + electron_mass_c2) /
                                   (deltaMomentum * totalMomentum);
  G4double sint = 1.0 - cost*cost;
  if(sint > 0.0) sint = std::sqrt(sint);


  G4ThreeVector deltaDirection(-sint*std::cos(phi),-sint*std::sin(phi), cost) ;
  deltaDirection.rotateUz(direction);

  // primary change
  kineticEnergy -= deltaKinEnergy;
  fParticleChange->SetProposedKineticEnergy(kineticEnergy);

  if(kineticEnergy > DBL_MIN) {
    G4ThreeVector dir = totalMomentum*direction - deltaMomentum*deltaDirection;
    direction = dir.unit();
    fParticleChange->SetProposedMomentumDirection(direction);
  }

  // create G4DynamicParticle object for delta ray
  G4DynamicParticle* delta = new G4DynamicParticle(theElectron,deltaDirection,deltaKinEnergy);
  vdp->push_back(delta);

  // get interaction frame
  G4ThreeVector  nInteractionFrame = 
    G4PolarizationHelper::GetFrame(direction,deltaDirection);

  if (crossSectionCalculator) {
    // calculate mean final state polarizations

    theBeamPolarization.InvRotateAz(nInteractionFrame,direction);
    theTargetPolarization.InvRotateAz(nInteractionFrame,direction);
    crossSectionCalculator->Initialize(x,gam,phi,theBeamPolarization,
				       theTargetPolarization,2);

    // electron/positron
    fPositronPolarization=crossSectionCalculator->GetPol2();
    fPositronPolarization.RotateAz(nInteractionFrame,direction);

    fParticleChange->ProposePolarization(fPositronPolarization);

    // electron
    fElectronPolarization=crossSectionCalculator->GetPol3();
    fElectronPolarization.RotateAz(nInteractionFrame,deltaDirection);
    delta->SetPolarization(fElectronPolarization.x(),
			   fElectronPolarization.y(),
			   fElectronPolarization.z());
  }
  else {
    fPositronPolarization=G4ThreeVector();
    fElectronPolarization=G4ThreeVector();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
