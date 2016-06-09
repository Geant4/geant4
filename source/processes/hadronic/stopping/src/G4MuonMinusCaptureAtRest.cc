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
// ------------------------------------------------------------
//      GEANT 4 class file
//
//      ------------ G4MuonMinusCaptureAtRest physics process ------
//                   by Larry Felawka (TRIUMF)
//                     E-mail: felawka@alph04.triumf.ca
//                   and Art Olin (TRIUMF)
//                     E-mail: olin@triumf.ca
//                            April 1998
//-----------------------------------------------------------------------------
//
// Modifications: 
// 18/08/2000  V.Ivanchenko Update description
// 17/05/2006  V.Ivanchenko Cleanup
//
//-----------------------------------------------------------------------------

#include "G4MuonMinusCaptureAtRest.hh"
#include "G4DynamicParticle.hh"
//#include "G4ParticleTypes.hh"
#include "Randomize.hh"
//#include <string.h>
#include <cmath>
//#include <stdio.h>
//#include <sys/types.h>
//#include <sys/stat.h>
#include "G4He3.hh"
#include "G4NeutrinoMu.hh"
#include "G4Fragment.hh"
#include "G4ReactionProductVector.hh"
#include "G4Proton.hh"
#include "G4PionPlus.hh"
#include "G4MuonMinus.hh"

G4MuonMinusCaptureAtRest::G4MuonMinusCaptureAtRest(const G4String& processName,
						   G4ProcessType   aType ) :
    G4VRestProcess (processName, aType), nCascade(0), targetZ(0),targetA(0)
{
  Cascade     = new G4GHEKinematicsVector [17];
  pSelector = new G4StopElementSelector();
  pEMCascade = new G4MuMinusCaptureCascade();
}

G4MuonMinusCaptureAtRest::~G4MuonMinusCaptureAtRest()
{
  delete [] Cascade;
  delete pSelector;
  delete pEMCascade;
}

G4bool G4MuonMinusCaptureAtRest::
IsApplicable(const G4ParticleDefinition& particle)
{
  return ( &particle == G4MuonMinus::MuonMinus() );

}

//
// Handles MuonMinuss at rest; a MuonMinus can either create secondaries or
// do nothing (in which case it should be sent back to decay-handling
// section
//
G4VParticleChange* G4MuonMinusCaptureAtRest::
AtRestDoIt(const G4Track& track,const G4Step&)
{
  aParticleChange.Initialize(track);

  // select element and get Z,A.
  G4Element* aEle = pSelector->GetElement(track.GetMaterial());
  targetZ = aEle->GetZ();
  targetA = aEle->GetN();

  G4IsotopeVector* isv = aEle->GetIsotopeVector();
  G4int ni = 0;
  if(isv) ni = isv->size();
  if(ni == 1) {
    targetA = G4double(aEle->GetIsotope(0)->GetN());
  } else if(ni > 0) {
    G4double* ab = aEle->GetRelativeAbundanceVector();
    G4double y = G4UniformRand();
    G4int j = -1;
    ni--;
    do {
      j++;
      y -= ab[j];
    } while (y > 0.0 && j < ni);
    targetA = G4double(aEle->GetIsotope(j)->GetN());
  }
  
  // Do the electromagnetic cascade of the muon in the nuclear field.
  nCascade = 0;
  targetMass = G4NucleiProperties::GetNuclearMass(targetA, targetZ);
  nCascade   = pEMCascade->DoCascade(targetZ, targetMass, Cascade);

  // Decide on Decay or Capture, and doit.
  G4double lambdac  = pSelector->GetMuonCaptureRate(targetZ, targetA);
  G4double lambdad  = pSelector->GetMuonDecayRate(targetZ, targetA);
  G4double lambda   = lambdac + lambdad;

  // ===  Throw for capture  time.

  G4double tDelay = -std::log(G4UniformRand()) / lambda;
  
  G4ReactionProductVector * captureResult = 0;
  G4int nEmSecondaries = nCascade;
  G4int nSecondaries = nCascade;

  if( G4UniformRand()*lambda > lambdac) 
  {
    pEMCascade->DoBoundMuonMinusDecay(targetZ, targetMass, &nEmSecondaries, Cascade);
  } 
  else 
  {
    captureResult = DoMuCapture();
  }

  // fill the final state
  if(captureResult) nSecondaries += captureResult->size();
  else nSecondaries = nEmSecondaries;
  aParticleChange.SetNumberOfSecondaries( nSecondaries );

  G4double globalTime = track.GetGlobalTime();
  G4ThreeVector position = track.GetPosition();
  // Store nuclear cascade
  if(captureResult)
  {
    for ( size_t isec = 0; isec < captureResult->size(); isec++ ) 
    {
      G4ReactionProduct* aParticle = captureResult->operator[](isec);
      G4DynamicParticle * aNewParticle = new G4DynamicParticle();
      aNewParticle->SetDefinition( aParticle->GetDefinition() );
      G4LorentzVector itV(aParticle->GetTotalEnergy(), aParticle->GetMomentum());
      aNewParticle->SetMomentum(itV.vect());
      G4double localtime = globalTime + tDelay + aParticle->GetTOF();
      G4Track* aNewTrack = new G4Track( aNewParticle, localtime, position);
	  	aNewTrack->SetTouchableHandle(track.GetTouchableHandle());
      aParticleChange.AddSecondary( aNewTrack );
    }
  }

  // Store electromagnetic cascade

  if(nEmSecondaries > 0) {

    for ( G4int isec = 0; isec < nEmSecondaries; isec++ ) {
      G4ParticleDefinition* pd = Cascade[isec].GetParticleDef();
      G4double localtime = globalTime;
      if(isec >= nCascade) localtime += tDelay;
      if(pd) {
        G4DynamicParticle* aNewParticle = new G4DynamicParticle;
        aNewParticle->SetDefinition( pd );
        aNewParticle->SetMomentum( Cascade[isec].GetMomentum() );

        G4Track* aNewTrack = new G4Track( aNewParticle, localtime, position );
		    aNewTrack->SetTouchableHandle(track.GetTouchableHandle());
        aParticleChange.AddSecondary( aNewTrack );
      }
    }
  }

  aParticleChange.ProposeLocalEnergyDeposit(0.0);
  aParticleChange.ProposeTrackStatus(fStopAndKill); 

  return &aParticleChange;
}

G4ReactionProductVector * G4MuonMinusCaptureAtRest::DoMuCapture()
{
  static G4double zeff[100] = {
    1.,1.98,2.95,3.89,4.8,5.72,6.61,7.49,8.32,9.12,9.95,10.69,11.48,12.22,
    12.91,13.64,14.24,14.89,15.53,16.15,16.75,17.38,18.04,18.49,
    19.06,19.59,20.1,20.66,21.12,21.61,22.02,22.43,22.84,23.24,
    23.65,24.06,24.47,24.85,25.23,25.61,25.99,26.37,26.69,27.,
    27.32,27.63,27.95,28.2,28.42,28.64,28.79,29.03,29.27,29.51,
    29.75,29.99,30.2,30.36,30.53,30.69,30.85,31.01,31.18,31.34,
    31.48,31.62,31.76,31.9,32.05,32.19,32.33,32.47,32.61,32.76,
    32.94,33.11,33.29,33.46,33.64,33.81,34.21,34.18,34.,34.1,
    34.21,34.31,34.42,34.52,34.63,34.73,34.84,34.94,35.04,35.15,
    35.25,35.36,35.46,35.57,35.67,35.78 };

  // Get the muon 4-vector
  //  G4cout << "G4MuonMinusCaptureAtRest::DoMuCapture called " << G4endl;
  G4int idxx = G4lrint(targetZ)-1;
  if(idxx>99) idxx=99;  
  G4double q = zeff[idxx];
  G4double zeff2 = q*q;
  G4double mumass = G4MuonMinus::MuonMinus()->GetPDGMass();
  G4double muonBindingEnergy = 0.5*zeff2*mumass*fine_structure_const*fine_structure_const;
  // Energy on K-shell
  G4double muEnergy = mumass + muonBindingEnergy;
  G4double availableEnergy = targetMass + mumass - muonBindingEnergy;

  G4double cost = 2.*G4UniformRand() - 1.0;
  G4double sint = std::sqrt((1.0 - cost)*(1.0 + cost));
  G4double phi  = twopi*G4UniformRand();
  G4ThreeVector aMu3Mom(sint*std::cos(phi),sint*std::sin(phi),cost);
  G4double pmu  = std::sqrt(muEnergy*muEnergy + mumass*mumass);
  G4LorentzVector aMuMom(muEnergy,aMu3Mom*pmu);

  G4double residualMass = G4NucleiProperties::GetNuclearMass(targetA, targetZ - 1.0);
  G4LorentzVector momResidual;
  G4ReactionProductVector * aPreResult;
  G4ReactionProduct* aNu = new G4ReactionProduct();
  aNu->SetDefinition( G4NeutrinoMu::NeutrinoMu() );

  // proton as a target
  if(targetA < 1.5) {

    G4double Ecms = mumass + proton_mass_c2 - muonBindingEnergy;
    G4double Enu  = 0.5*(Ecms - neutron_mass_c2*neutron_mass_c2/Ecms);

    // make the nu, and transform to lab;
    G4double cost = 2.*G4UniformRand() - 1.0;
    G4double sint = std::sqrt((1.0 - cost)*(1.0 + cost));
    G4double phi  = twopi*G4UniformRand();
    G4ThreeVector nu3Mom(sint*std::cos(phi),sint*std::sin(phi),cost);
    nu3Mom *= Enu;

    aPreResult = new G4ReactionProductVector();

    G4ReactionProduct* aN = new G4ReactionProduct();
    aN->SetDefinition( G4Neutron::Neutron() );
    aN->SetTotalEnergy( Ecms - Enu );
    aN->SetMomentum( -nu3Mom );

    aNu->SetTotalEnergy( Enu );
    aNu->SetMomentum( nu3Mom );
    aPreResult->push_back(aN ); 
    aPreResult->push_back(aNu); 
    if(verboseLevel > 1)
      G4cout << "G4MuonMinusCaptureAtRest::DoMuCapture on H " 
	     <<" EkinN(MeV)= " << (Ecms - Enu - neutron_mass_c2)/GeV
	     <<" Enu(MeV)= "<<aNu->GetTotalEnergy()/MeV<<G4endl;

    return aPreResult;
  }

  // pick random proton inside nucleus 
  G4double eEx=0;
  do {
    theN.Init(targetA, targetZ); 
    G4ThreeVector fermiMom;
    G4Nucleon * aNucleon = 0;
    G4int theProtonCounter = G4lrint( 0.5 + targetZ * G4UniformRand() );
    G4int counter = 0;
    theN.StartLoop();

    while( (aNucleon=theN.GetNextNucleon()) ) {

      if( aNucleon->GetDefinition() == G4Proton::Proton() ) {
	counter++;
	if(counter == theProtonCounter) {
          fermiMom = aNucleon->GetMomentum().vect();
	  break;
	}
      }
    }

    // Get the nu momentum in the CMS
    G4LorentzVector thePMom(std::sqrt(proton_mass_c2*proton_mass_c2 + fermiMom.mag2()), 
			    fermiMom);
    G4LorentzVector   theCMS = thePMom + aMuMom;
    G4ThreeVector bst = theCMS.boostVector();

    momResidual = G4LorentzVector(0.0,0.0,0.0,availableEnergy);

    G4double Ecms = theCMS.mag();
    G4double Enu  = 0.5*(Ecms - neutron_mass_c2*neutron_mass_c2/Ecms);

    // make the nu, and transform to lab;
    cost = 2.*G4UniformRand() - 1.0;
    sint = std::sqrt((1.0 - cost)*(1.0 + cost));
    phi  = twopi*G4UniformRand();
    G4ThreeVector nu3Mom(sint*std::cos(phi),sint*std::sin(phi),cost);
    G4LorentzVector nuMom(Enu,aMu3Mom*Enu);
    
    // make the neutrino an mu-neutrino with the above momentum and get the residual properties
    nuMom.boost(bst);
    momResidual -= nuMom;

    // nu in lab.
    aNu->SetTotalEnergy( nuMom.t() );
    aNu->SetMomentum( nuMom.vect() );

    // Call pre-compound on the rest.
    eEx = momResidual.mag() - residualMass;
  } while(eEx <= 0.0);
  
  G4ThreeVector fromBreit = momResidual.boostVector();
  G4LorentzVector fscm(0.0,0.0,0.0,  momResidual.mag());
  G4Fragment anInitialState;
  anInitialState.SetA(G4lrint(targetA));
  anInitialState.SetZ(G4lrint(targetZ) - 1);
  anInitialState.SetNumberOfParticles(2);
  anInitialState.SetNumberOfCharged(0);
  anInitialState.SetNumberOfHoles(1);
  anInitialState.SetMomentum(fscm);
  aPreResult = theHandler.BreakItUp(anInitialState);

  G4ReactionProductVector::iterator ires;
  G4double eBal = availableEnergy;
  for(ires=aPreResult->begin(); ires!=aPreResult->end(); ires++)
  {
    G4LorentzVector itV((*ires)->GetTotalEnergy(), (*ires)->GetMomentum());
    itV.boost(fromBreit);
    (*ires)->SetTotalEnergy(itV.t());
    (*ires)->SetMomentum(itV.vect());
    eBal -= itV.t();
  }
  // fill neutrino into result
  aPreResult->push_back(aNu);
  eBal -= aNu->GetTotalEnergy();
  if(verboseLevel > 1)
    G4cout << "G4MuonMinusCaptureAtRest::DoMuCapture:  Nsec= " 
	   << aPreResult->size() << " Ebalance(MeV)= " << eBal/MeV
           << " Eex(MeV)= " << eEx/MeV
	   <<" E0(GeV)= " <<availableEnergy/GeV
	   <<" Enu(MeV)= "<<aNu->GetTotalEnergy()/MeV<<G4endl;

  return aPreResult;
} 

