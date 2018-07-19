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
//
// $Id: G4PolarizedComptonModel.cc 96114 2016-03-16 18:51:33Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4PolarizedComptonModel
//
// Author:        Andreas Schaelicke
//
// Creation date: 01.05.2005
//
// Modifications:
// 18-07-06 use newly calculated cross sections (P. Starovoitov)
// 21-08-05 update interface (A. Schaelicke)
//
// Class Description:
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4PolarizedComptonModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "Randomize.hh"
#include "G4DataVector.hh"
#include "G4ParticleChangeForGamma.hh"

#include "G4StokesVector.hh"
#include "G4PolarizationManager.hh"
#include "G4PolarizationHelper.hh"
#include "G4PolarizedComptonCrossSection.hh"

#include "G4SystemOfUnits.hh"
#include "G4Log.hh"
#include "G4Exp.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

static const G4int nlooplim = 10000;

G4PolarizedComptonModel::G4PolarizedComptonModel(const G4ParticleDefinition*,
						 const G4String& nam)
  : G4KleinNishinaCompton(nullptr,nam),
    verboseLevel(0)
{
  crossSectionCalculator = new G4PolarizedComptonCrossSection();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PolarizedComptonModel::~G4PolarizedComptonModel()
{
  delete crossSectionCalculator;
}

G4double G4PolarizedComptonModel::ComputeAsymmetryPerAtom
                       (G4double gammaEnergy, G4double /*Z*/)
 
{
  G4double asymmetry = 0.0 ;

  G4double k0 = gammaEnergy / electron_mass_c2 ;
  G4double k1 = 1. + 2.*k0 ;

  asymmetry = -k0;
  asymmetry *= (k0 + 1.)*sqr(k1)*G4Log(k1) - 2.*k0*(5.*sqr(k0) + 4.*k0 + 1.);
  asymmetry /= ((k0 - 2.)*k0  -2.)*sqr(k1)*G4Log(k1) + 2.*k0*(k0*(k0 + 1.)*(k0 + 8.) + 2.);		

  // G4cout<<"energy = "<<GammaEnergy<<"  asymmetry = "<<asymmetry<<"\t\t GAM = "<<k0<<G4endl;
  if (asymmetry>1.) G4cout<<"ERROR in G4PolarizedComptonModel::ComputeAsymmetryPerAtom"<<G4endl;

  return asymmetry;
}


G4double G4PolarizedComptonModel::ComputeCrossSectionPerAtom(
                                const G4ParticleDefinition* pd,
                                      G4double kinEnergy, 
                                      G4double Z, 
                                      G4double A, 
                                      G4double cut,
                                      G4double emax)
{
  double xs = 
    G4KleinNishinaCompton::ComputeCrossSectionPerAtom(pd,kinEnergy,
						      Z,A,cut,emax);
  G4double polzz = theBeamPolarization.p3()*theTargetPolarization.z();
  if (polzz > 0.0) {
    G4double asym = ComputeAsymmetryPerAtom(kinEnergy, Z);  
    xs *= (1.+polzz*asym);
  }
  return xs;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PolarizedComptonModel::SampleSecondaries(
                              std::vector<G4DynamicParticle*>* fvect,
                              const G4MaterialCutsCouple*,
			      const G4DynamicParticle* aDynamicGamma,
			      G4double, G4double)
{
  // do nothing below the threshold
  if(aDynamicGamma->GetKineticEnergy() <= LowEnergyLimit()) { return; }

  const G4Track * aTrack = fParticleChange->GetCurrentTrack();
  G4VPhysicalVolume*  aPVolume  = aTrack->GetVolume();
  G4LogicalVolume*    aLVolume  = aPVolume->GetLogicalVolume();

  if (verboseLevel >= 1) {
    G4cout<<"G4PolarizedComptonModel::SampleSecondaries in "
          <<  aLVolume->GetName() <<G4endl;
  }
  G4PolarizationManager * polarizationManager = 
    G4PolarizationManager::GetInstance();

  // obtain polarization of the beam
  theBeamPolarization =  aDynamicGamma->GetPolarization();
  theBeamPolarization.SetPhoton();

  // obtain polarization of the media
  G4bool targetIsPolarized = polarizationManager->IsPolarized(aLVolume);
  theTargetPolarization = 
    polarizationManager->GetVolumePolarization(aLVolume);

  // if beam is linear polarized or target is transversely polarized 
  // determine the angle to x-axis
  // (assumes same PRF as in the polarization definition)

  G4ThreeVector gamDirection0 = aDynamicGamma->GetMomentumDirection();

  // transfere theTargetPolarization 
  // into the gamma frame (problem electron is at rest)
  if (targetIsPolarized) {
    theTargetPolarization.rotateUz(gamDirection0);
  }
  // The scattered gamma energy is sampled according to 
  // Klein - Nishina formula.
  // The random number techniques of Butcher & Messel are used 
  // (Nuc Phys 20(1960),15).
  // Note : Effects due to binding of atomic electrons are negliged.
 
  G4double gamEnergy0 = aDynamicGamma->GetKineticEnergy();
  G4double E0_m = gamEnergy0 / electron_mass_c2 ;

  //
  // sample the energy rate of the scattered gamma 
  //

  G4double epsilon, sint2;
  G4double onecost = 0.0;
  G4double Phi     = 0.0;
  G4double greject = 1.0;
  G4double cosTeta = 1.0;
  G4double sinTeta = 0.0;

  G4double eps0       = 1./(1. + 2.*E0_m);
  G4double epsilon0sq = eps0*eps0;
  G4double alpha1     = - G4Log(eps0);
  G4double alpha2     = alpha1 + 0.5*(1.- epsilon0sq);

  G4double polarization = 
    theBeamPolarization.p3()*theTargetPolarization.p3();

  CLHEP::HepRandomEngine* rndmEngineMod = G4Random::getTheEngine();
  G4int nloop = 0;
  G4bool end = false;

  G4double rndm[3];

  do {
    do {
      ++nloop;
      // false interaction if too many iterations
      if(nloop > nlooplim) { 
	PrintWarning(aDynamicGamma, nloop, greject, onecost, Phi, 
		     "too many iterations"); 
	return; 
      }

      // 3 random numbers to sample scattering
      rndmEngineMod->flatArray(3, rndm);

      if ( alpha1 > alpha2*rndm[0]) {
	epsilon   = G4Exp(-alpha1*rndm[1]);   // epsilon0**r
      } else {
	epsilon = std::sqrt(epsilon0sq + (1.- epsilon0sq)*rndm[1]);
      }

      onecost = (1.- epsilon)/(epsilon*E0_m);
      sint2   = onecost*(2.-onecost);

      G4double gdiced = 2.*(1./epsilon+epsilon);
      G4double gdist  = 1./epsilon + epsilon - sint2 
	- polarization*(1./epsilon-epsilon)*(1.-onecost);

      greject = gdist/gdiced;

      if (greject > 1.0) { 
	PrintWarning(aDynamicGamma, nloop, greject, onecost, Phi, 
		     "theta majoranta wrong"); 
      }
      // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
    } while (greject < rndm[2]);

    // assuming phi loop sucessful
    end = true;
 
    //
    // scattered gamma angles. ( Z - axis along the parent gamma)
    //
    cosTeta = 1. - onecost; 
    sinTeta = std::sqrt(sint2);
    do {
      ++nloop;

      // 2 random numbers to sample scattering
      rndmEngineMod->flatArray(2, rndm);

      // false interaction if too many iterations
      Phi = twopi * rndm[0];
      if(nloop > nlooplim) { 
	PrintWarning(aDynamicGamma, nloop, greject, onecost, Phi, 
		     "too many iterations"); 
	return; 
      }

      G4double gdiced = 1./epsilon + epsilon - sint2 
	+ std::abs(theBeamPolarization.p3())*
	( std::abs((1./epsilon-epsilon)*cosTeta*theTargetPolarization.p3())
	  +(1.-epsilon)*sinTeta*(std::sqrt(sqr(theTargetPolarization.p1()) 
					   + sqr(theTargetPolarization.p2()))))
	+sint2*(std::sqrt(sqr(theBeamPolarization.p1()) + 
			  sqr(theBeamPolarization.p2())));

      G4double gdist = 1./epsilon + epsilon - sint2 
	+ theBeamPolarization.p3()*
	((1./epsilon-epsilon)*cosTeta*theTargetPolarization.p3()
	 +(1.-epsilon)*sinTeta*(std::cos(Phi)*theTargetPolarization.p1()+
				std::sin(Phi)*theTargetPolarization.p2()))
	-sint2*(std::cos(2.*Phi)*theBeamPolarization.p1()
		+std::sin(2.*Phi)*theBeamPolarization.p2());
      greject = gdist/gdiced;

      if (greject > 1.0) {
	PrintWarning(aDynamicGamma, nloop, greject, onecost, Phi, 
		     "phi majoranta wrong"); 
      }

      if(greject < 1.e-3) {
	PrintWarning(aDynamicGamma, nloop, greject, onecost, Phi, 
		     "phi loop ineffective"); 
	// restart theta loop
        end = false;
        break;
      }
     
      // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
    } while (greject < rndm[1]);
  } while(!end);
  G4double dirx = sinTeta*std::cos(Phi), diry = sinTeta*std::sin(Phi), 
    dirz = cosTeta;

  //
  // update G4VParticleChange for the scattered gamma
  //
   
  G4ThreeVector gamDirection1 ( dirx,diry,dirz );
  gamDirection1.rotateUz(gamDirection0);
  G4double gamEnergy1 = epsilon*gamEnergy0;

  G4double edep = 0.0;
  if(gamEnergy1 > lowestSecondaryEnergy) {
    fParticleChange->ProposeMomentumDirection(gamDirection1);
    fParticleChange->SetProposedKineticEnergy(gamEnergy1);
  } else { 
    fParticleChange->ProposeTrackStatus(fStopAndKill);
    fParticleChange->SetProposedKineticEnergy(0.0);
    edep = gamEnergy1;
  }
 
  // 
  // calculate Stokesvector of final state photon and electron
  //
  G4ThreeVector  nInteractionFrame = 
    G4PolarizationHelper::GetFrame(gamDirection1,gamDirection0);

  // transfere theBeamPolarization and theTargetPolarization 
  // into the interaction frame (note electron is in gamma frame)
  if (verboseLevel>=1) {
    G4cout << "========================================\n";
    G4cout << " nInteractionFrame = " <<nInteractionFrame<<"\n";
    G4cout << " GammaDirection0 = " <<gamDirection0<<"\n";
    G4cout << " gammaPolarization = " <<theBeamPolarization<<"\n";
    G4cout << " electronPolarization = " <<theTargetPolarization<<"\n";
  }

  theBeamPolarization.InvRotateAz(nInteractionFrame,gamDirection0);
  theTargetPolarization.InvRotateAz(nInteractionFrame,gamDirection0);

  if (verboseLevel>=1) {
    G4cout << "----------------------------------------\n";
    G4cout << " gammaPolarization = " <<theBeamPolarization<<"\n";
    G4cout << " electronPolarization = " <<theTargetPolarization<<"\n";
    G4cout << "----------------------------------------\n";
  }

  // initialize the polarization transfer matrix
  crossSectionCalculator->Initialize(epsilon,E0_m,0.,
				     theBeamPolarization,
				     theTargetPolarization,2);
  
  if(gamEnergy1 > lowestSecondaryEnergy) {
 
    // in interaction frame
    // calculate polarization transfer to the photon (in interaction plane)
    finalGammaPolarization = crossSectionCalculator->GetPol2();
    if (verboseLevel>=1) {
      G4cout << " gammaPolarization1 = " <<finalGammaPolarization<<"\n";
    }
    finalGammaPolarization.SetPhoton();

    // translate polarization into particle reference frame
    finalGammaPolarization.RotateAz(nInteractionFrame,gamDirection1);
    if (finalGammaPolarization.mag() > 1.+1.e-8){
      G4cout<<"ERROR in Polarizaed Compton Scattering !"<<G4endl;
      G4cout<<"Polarization of final photon more than 100%"<<G4endl;
      G4cout<<finalGammaPolarization<<" mag = "
	    <<finalGammaPolarization.mag()<<G4endl;
    }
    //store polarization vector
    fParticleChange->ProposePolarization(finalGammaPolarization);
    if (verboseLevel>=1) {
      G4cout << " gammaPolarization1 = " <<finalGammaPolarization<<"\n";
      G4cout << " GammaDirection1 = " <<gamDirection1<<"\n";
    }
  }

  //
  // kinematic of the scattered electron
  //
  G4double eKinEnergy = gamEnergy0 - gamEnergy1;

  if (eKinEnergy > lowestSecondaryEnergy) {
  
    G4ThreeVector eDirection = 
      gamEnergy0*gamDirection0 - gamEnergy1*gamDirection1;
    eDirection = eDirection.unit();

    finalElectronPolarization = crossSectionCalculator->GetPol3();
    if (verboseLevel>=1) {
      G4cout << " electronPolarization1 = " 
	     <<finalElectronPolarization<<"\n";
    }
    // transfer into particle reference frame
    finalElectronPolarization.RotateAz(nInteractionFrame,eDirection);
    if (verboseLevel>=1) {
      G4cout << " electronPolarization1 = " 
	     <<finalElectronPolarization<<"\n";
      G4cout << " ElecDirection = " <<eDirection<<"\n";
    }

    // create G4DynamicParticle object for the electron.
    G4DynamicParticle* aElectron = 
      new G4DynamicParticle(theElectron,eDirection,eKinEnergy);
    //store polarization vector
    if (finalElectronPolarization.mag() > 1.+1.e-8){
      G4cout<<"ERROR in Polarizaed Compton Scattering !"<<G4endl;
      G4cout<<"Polarization of final electron more than 100%"<<G4endl;
      G4cout<<finalElectronPolarization<<" mag = "
	    <<finalElectronPolarization.mag()<<G4endl;
    }
    aElectron->SetPolarization(finalElectronPolarization.p1(),
			       finalElectronPolarization.p2(),
			       finalElectronPolarization.p3());
    fvect->push_back(aElectron);
  } else {
    edep += eKinEnergy;  
  }
  // energy balance
  if(edep > 0.0) { 
    fParticleChange->ProposeLocalEnergyDeposit(edep);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void 
G4PolarizedComptonModel::PrintWarning(const G4DynamicParticle* dp, G4int nloop,
				      G4double grej, G4double onecos, 
				      G4double phi, const G4String sss) const
{
  G4ExceptionDescription ed;
  ed << "Problem of scattering sampling: " << sss << "\n"
     << "Niter= " << nloop << " grej= " << grej << " cos(theta)= " 
     << 1.0-onecos << " phi= " << phi << "\n"
     << "Gamma E(MeV)= " << dp->GetKineticEnergy()/MeV
     << " dir= " << dp->GetMomentumDirection() 
     << " pol= " << dp->GetPolarization();
  G4Exception("G4PolarizedComptonModel::SampleSecondaries","em0044",
	      JustWarning, ed, "");
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


