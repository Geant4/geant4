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
// $Id: G4PolarizedComptonModel.cc,v 1.4 2007-05-23 08:52:20 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "Randomize.hh"
#include "G4DataVector.hh"
#include "G4ParticleChangeForGamma.hh"


#include "G4StokesVector.hh"
#include "G4PolarizationManager.hh"
#include "G4PolarizationHelper.hh"
#include "G4PolarizedComptonCrossSection.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PolarizedComptonModel::G4PolarizedComptonModel(const G4ParticleDefinition*,
                                             const G4String& nam)
  : G4KleinNishinaCompton(0,nam),
    verboseLevel(0)
{
  crossSectionCalculator=new G4PolarizedComptonCrossSection();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PolarizedComptonModel::~G4PolarizedComptonModel()
{
  if (crossSectionCalculator) delete crossSectionCalculator;
}



G4double G4PolarizedComptonModel::ComputeAsymmetryPerAtom
                       (G4double gammaEnergy, G4double /*Z*/)
 
{
 G4double asymmetry = 0.0 ;

 G4double k0 = gammaEnergy / electron_mass_c2 ;
 G4double k1 = 1 + 2*k0 ;

 asymmetry = -k0;
 asymmetry *= (k0 + 1.)*sqr(k1)*std::log(k1) - 2.*k0*(5.*sqr(k0) + 4.*k0 + 1.);
 asymmetry /= ((k0 - 2.)*k0  -2.)*sqr(k1)*std::log(k1) + 2.*k0*(k0*(k0 + 1.)*(k0 + 8.) + 2.);		

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
  if (polzz!=0) {
    G4double asym=ComputeAsymmetryPerAtom(kinEnergy, Z);  
    xs*=(1.+polzz*asym);
  }
  return xs;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PolarizedComptonModel::SampleSecondaries(std::vector<G4DynamicParticle*>* fvect,
						const G4MaterialCutsCouple*,
						const G4DynamicParticle* aDynamicGamma,
						G4double,
						G4double)
{
  const G4Track * aTrack = fParticleChange->GetCurrentTrack();
  G4VPhysicalVolume*  aPVolume  = aTrack->GetVolume();
  G4LogicalVolume*    aLVolume  = aPVolume->GetLogicalVolume();

  if (verboseLevel>=1) 
    G4cout<<"G4PolarizedComptonModel::SampleSecondaries in "
          <<  aLVolume->GetName() <<G4endl;

  G4PolarizationManager * polarizationManager = G4PolarizationManager::GetInstance();

  // obtain polarization of the beam
  theBeamPolarization =  aDynamicGamma->GetPolarization();
  theBeamPolarization.SetPhoton();

  // obtain polarization of the media
  const G4bool targetIsPolarized = polarizationManager->IsPolarized(aLVolume);
  theTargetPolarization = polarizationManager->GetVolumePolarization(aLVolume);

  // if beam is linear polarized or target is transversely polarized 
  // determine the angle to x-axis
  // (assumes same PRF as in the polarization definition)

  G4ThreeVector gamDirection0 = aDynamicGamma->GetMomentumDirection();

  // transfere theTargetPolarization 
  // into the gamma frame (problem electron is at rest)
  if (targetIsPolarized)
    theTargetPolarization.rotateUz(gamDirection0);

  // The scattered gamma energy is sampled according to Klein - Nishina formula.
  // The random number techniques of Butcher & Messel are used 
  // (Nuc Phys 20(1960),15).
  // Note : Effects due to binding of atomic electrons are negliged.
 
  G4double gamEnergy0 = aDynamicGamma->GetKineticEnergy();
  G4double E0_m = gamEnergy0 / electron_mass_c2 ;


  //
  // sample the energy rate of the scattered gamma 
  //

  G4double epsilon, epsilonsq, onecost, sint2, greject ;

  G4double epsilon0   = 1./(1. + 2.*E0_m);
  G4double epsilon0sq = epsilon0*epsilon0;
  G4double alpha1     = - std::log(epsilon0);
  G4double alpha2     = 0.5*(1.- epsilon0sq);

  G4double polarization = theBeamPolarization.p3()*theTargetPolarization.p3();
  do {
    if ( alpha1/(alpha1+alpha2) > G4UniformRand() ) {
      epsilon   = std::exp(-alpha1*G4UniformRand());   // epsilon0**r
      epsilonsq = epsilon*epsilon; 

    } else {
      epsilonsq = epsilon0sq + (1.- epsilon0sq)*G4UniformRand();
      epsilon   = std::sqrt(epsilonsq);
    };

    onecost = (1.- epsilon)/(epsilon*E0_m);
    sint2   = onecost*(2.-onecost);


    G4double gdiced = 2.*(1./epsilon+epsilon);
    G4double gdist  = 1./epsilon + epsilon - sint2 
      - polarization*(1./epsilon-epsilon)*(1.-onecost);

    greject = gdist/gdiced;

    if (greject>1) G4cout<<"ERROR in PolarizedComptonScattering::PostStepDoIt\n"
			 <<" costh rejection does not work properly: "<<greject<<G4endl;

  } while (greject < G4UniformRand());
 
  //
  // scattered gamma angles. ( Z - axis along the parent gamma)
  //

  G4double cosTeta = 1. - onecost; 
  G4double sinTeta = std::sqrt (sint2);
  G4double Phi;
  do {
    Phi     = twopi * G4UniformRand();
     G4double gdiced = 1./epsilon + epsilon - sint2 
       + std::abs(theBeamPolarization.p3())*
       ( std::abs((1./epsilon-epsilon)*cosTeta*theTargetPolarization.p3())
	+(1.-epsilon)*sinTeta*(std::sqrt(sqr(theTargetPolarization.p1()) 
				    + sqr(theTargetPolarization.p2()))))
       +sint2*(std::sqrt(sqr(theBeamPolarization.p1()) + sqr(theBeamPolarization.p2())));

     G4double gdist = 1./epsilon + epsilon - sint2 
       + theBeamPolarization.p3()*
       ((1./epsilon-epsilon)*cosTeta*theTargetPolarization.p3()
	+(1.-epsilon)*sinTeta*(std::cos(Phi)*theTargetPolarization.p1()+
			       std::sin(Phi)*theTargetPolarization.p2()))
       -sint2*(std::cos(2.*Phi)*theBeamPolarization.p1()
	       +std::sin(2.*Phi)*theBeamPolarization.p2());
     greject = gdist/gdiced;

    if (greject>1.+1.e-10 || greject<0) G4cout<<"ERROR in PolarizedComptonScattering::PostStepDoIt\n"
				      <<" phi rejection does not work properly: "<<greject<<G4endl;

    if (greject<1.e-3) {
      G4cout<<"ERROR in PolarizedComptonScattering::PostStepDoIt\n"
	    <<" phi rejection does not work properly: "<<greject<<"\n";
      G4cout<<" greject="<<greject<<"  phi="<<Phi<<"   cost="<<cosTeta<<"\n";
      G4cout<<" gdiced="<<gdiced<<"   gdist="<<gdist<<"\n";
      G4cout<<" eps="<<epsilon<<"    1/eps="<<1./epsilon<<"\n";
    }
     
  } while (greject < G4UniformRand());
  G4double dirx = sinTeta*std::cos(Phi), diry = sinTeta*std::sin(Phi), dirz = cosTeta;

  //
  // update G4VParticleChange for the scattered gamma
  //
   
  G4ThreeVector gamDirection1 ( dirx,diry,dirz );
  gamDirection1.rotateUz(gamDirection0);
  G4double gamEnergy1 = epsilon*gamEnergy0;
  fParticleChange->SetProposedKineticEnergy(gamEnergy1);


  if(gamEnergy1 > lowestGammaEnergy) {
    fParticleChange->ProposeMomentumDirection(gamDirection1);
  } else { 
    fParticleChange->ProposeTrackStatus(fStopAndKill);
    gamEnergy1 += fParticleChange->GetLocalEnergyDeposit();
    fParticleChange->ProposeLocalEnergyDeposit(gamEnergy1);
  }
 
  //
  // kinematic of the scattered electron
  //

  G4double eKinEnergy = gamEnergy0 - gamEnergy1;
  G4ThreeVector eDirection = gamEnergy0*gamDirection0 - gamEnergy1*gamDirection1;
  eDirection = eDirection.unit();

  // 
  // calculate Stokesvector of final state photon and electron
  //
  G4ThreeVector  nInteractionFrame;
  if((gamEnergy1 > lowestGammaEnergy) ||
     (eKinEnergy > DBL_MIN)) {

    // determine interaction plane
//     nInteractionFrame = 
//       G4PolarizationHelper::GetFrame(gamDirection1,eDirection);
    if (gamEnergy1 > lowestGammaEnergy) 
      nInteractionFrame = G4PolarizationHelper::GetFrame(gamDirection1,gamDirection0);
    else 
      nInteractionFrame = G4PolarizationHelper::GetFrame(gamDirection0, eDirection);

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
  }

  //  if(eKinEnergy > DBL_MIN)
  {
    // in interaction frame
    // calculate polarization transfer to the photon (in interaction plane)
    finalGammaPolarization = crossSectionCalculator->GetPol2();
    if (verboseLevel>=1) G4cout << " gammaPolarization1 = " <<finalGammaPolarization<<"\n";
    finalGammaPolarization.SetPhoton();

    // translate polarization into particle reference frame
    finalGammaPolarization.RotateAz(nInteractionFrame,gamDirection1);
    //store polarization vector
    fParticleChange->ProposePolarization(finalGammaPolarization);
    if (finalGammaPolarization.mag() > 1.+1.e-8){
      G4cout<<"ERROR in Polarizaed Compton Scattering !"<<G4endl;
      G4cout<<"Polarization of final photon more than 100%"<<G4endl;
      G4cout<<finalGammaPolarization<<" mag = "<<finalGammaPolarization.mag()<<G4endl;
    }
    if (verboseLevel>=1) {
      G4cout << " gammaPolarization1 = " <<finalGammaPolarization<<"\n";
      G4cout << " GammaDirection1 = " <<gamDirection1<<"\n";
    }
  }

  //    if (ElecKineEnergy > fminimalEnergy) {
  {
    finalElectronPolarization = crossSectionCalculator->GetPol3();
    if (verboseLevel>=1) 
      G4cout << " electronPolarization1 = " <<finalElectronPolarization<<"\n";

    // transfer into particle reference frame
    finalElectronPolarization.RotateAz(nInteractionFrame,eDirection);
    if (verboseLevel>=1) {
      G4cout << " electronPolarization1 = " <<finalElectronPolarization<<"\n";
      G4cout << " ElecDirection = " <<eDirection<<"\n";
    }
  }
  if (verboseLevel>=1)
    G4cout << "========================================\n";
       

  if(eKinEnergy > DBL_MIN) {

    // create G4DynamicParticle object for the electron.
    G4DynamicParticle* aElectron = new G4DynamicParticle(theElectron,eDirection,eKinEnergy);
    //store polarization vector
    if (finalElectronPolarization.mag() > 1.+1.e-8){
      G4cout<<"ERROR in Polarizaed Compton Scattering !"<<G4endl;
      G4cout<<"Polarization of final electron more than 100%"<<G4endl;
      G4cout<<finalElectronPolarization<<" mag = "<<finalElectronPolarization.mag()<<G4endl;
    }
    aElectron->SetPolarization(finalElectronPolarization.p1(),
			       finalElectronPolarization.p2(),
			       finalElectronPolarization.p3());
    fvect->push_back(aElectron);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


