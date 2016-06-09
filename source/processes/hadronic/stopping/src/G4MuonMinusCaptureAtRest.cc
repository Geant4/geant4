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
//

#include "G4MuonMinusCaptureAtRest.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleTypes.hh"
#include "Randomize.hh"
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "G4He3.hh"
#include "G4NeutrinoMu.hh"
#include "G4Fragment.hh"
#include "G4ReactionProductVector.hh"
#include "G4Proton.hh"
#include "G4PionPlus.hh"

G4MuonMinusCaptureAtRest::G4MuonMinusCaptureAtRest(const G4String& processName)
  : G4VRestProcess (processName), nCascade(0), targetCharge(0),
    targetAtomicMass(0)
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
/*
G4double G4MuonMinusCaptureAtRest::
GetMeanLifeTime(const G4Track& track, G4ForceCondition* )
{
  GetCaptureIsotope( track );
  return (tDelay );

}

void G4MuonMinusCaptureAtRest::
GetCaptureIsotope(const G4Track& track)
{

  // Ask selector to choose the element

  G4Material * aMaterial = track.GetMaterial();
  G4Element* theElement  = pSelector->GetElement(aMaterial);
  targetCharge           = theElement->GetZ();
  targetAtomicMass       = theElement->GetN();

  // Calculate total capture velosity

  G4double lambdac  = pSelector->GetMuonCaptureRate(targetCharge,targetAtomicMass);
           lambdac += pSelector->GetMuonDecayRate(targetCharge,targetAtomicMass);

  // ===  Throw for capture  time.

  G4double tDelay = -(G4double)log(G4UniformRand()) / lambdac;

  return;

} 
*/

G4double G4MuonMinusCaptureAtRest::
AtRestGetPhysicalInteractionLength(const G4Track&, G4ForceCondition* condition)
{
  // beggining of tracking
  //ResetNumberOfInteractionLengthLeft();

  // condition is set to "Not Forced"
  *condition = NotForced;
  // condition is set to "ExclusivelyForced" by V.Ivanchenko
  //     *condition = ExclusivelyForced;

  // get mean life time
  //  currentInteractionLength = GetMeanLifeTime(track, condition);
  /*
  if ((currentInteractionLength <0.0) || (verboseLevel>2)){
    G4cout << "G4MuonMinusCaptureAtRestProcess::AtRestGetPhysicalInteractionLength ";
    G4cout << "[ " << GetProcessName() << "]" <<G4endl;
    track.GetDynamicParticle()->DumpInfo();
    G4cout << " in Material  " << track.GetMaterial()->GetName() <<G4endl;
    G4cout << "MeanLifeTime = " << currentInteractionLength/ns << "[ns]" <<G4endl;
  }

  // Return 0 interaction length to get for this process
  // the 100% probability
  */
  return 0.0;
  //  return theNumberOfInteractionLengthLeft * currentInteractionLength;

}

//
// Handles MuonMinuss at rest; a MuonMinus can either create secondaries or
// do nothing (in which case it should be sent back to decay-handling
// section
//
G4VParticleChange* G4MuonMinusCaptureAtRest::
AtRestDoIt(const G4Track& track,const G4Step&)
{
  static int reactionCount=0;
  if(getenv("MuonMinusCaptureAtRestPrint"))
  {
    G4cout << " ++ G4MuonMinusCaptureAtRest::AtRestDoIt is start " << ++reactionCount<<G4endl;
  }
  aParticleChange.Initialize(track);

 // if (verboseLevel > 1) {
 // }

  // select element and get Z,A.
  G4Element* aEle = pSelector->GetElement(track.GetMaterial());
  targetCharge = aEle->GetZ();
  targetAtomicMass = aEle->GetN();
  
  // Do the electromagnetic cascade of the muon in the nuclear field.
  nCascade = 0;
  G4double mass = GetTargetMass(targetAtomicMass,targetCharge);
  nCascade      = pEMCascade->DoCascade(targetCharge, mass, Cascade);

  // Decide on Decay or Capture, and doit.
  G4double lambdac  = pSelector->GetMuonCaptureRate(targetCharge,targetAtomicMass);
  G4double lambdad  = pSelector->GetMuonDecayRate(targetCharge,targetAtomicMass);
  G4double lambda   = lambdac + lambdad;

  // ===  Throw for capture  time.

  G4double tDelay = -std::log(G4UniformRand()) / lambda;
  
  G4ReactionProductVector * captureResult=0;
  G4int nEmSecondaries = nCascade;
  G4int nSecondaries = nCascade;

  if( G4UniformRand()*lambda > lambdac) 
  {
    pEMCascade->DoBoundMuonMinusDecay(targetCharge, mass, &nEmSecondaries, Cascade);
  } 
  else 
  {
    captureResult = DoMuCapture(track.GetKineticEnergy());
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
      if(isec > nCascade) localtime += tDelay;
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

  aParticleChange.ProposeTrackStatus(fStopAndKill); // Kill the incident MuonMinus

//   clear InteractionLengthLeft

  ResetNumberOfInteractionLengthLeft();

  if(getenv("MuonMinusCaptureAtRestPrint"))
  {
    G4cout << " -- G4MuonMinusCaptureAtRest::AtRestDoIt is end " << reactionCount<<G4endl;
  }
  return &aParticleChange;

}


G4double G4MuonMinusCaptureAtRest::GetTargetMass(G4double a, G4double z)
{
  G4double result;
  result = G4ParticleTable::GetParticleTable()
           ->GetIonTable()->GetIonMass(G4lrint(z), G4lrint(a));
  return result;
} 



G4ReactionProductVector * G4MuonMinusCaptureAtRest::DoMuCapture(G4double aMuKinetic)
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
  G4int idxx = G4lrint(targetCharge)-1;
  if(idxx>99) idxx=99;  
  G4double q = zeff[idxx];
  G4double zeff2 = q*q;
  G4double muonBindingEnergy = 0.5*zeff2 * G4MuonMinus::MuonMinusDefinition()->GetPDGMass() 
                               * fine_structure_const*fine_structure_const;
  G4double availableEnergy = aMuKinetic;
  availableEnergy += G4MuonMinus::MuonMinusDefinition()->GetPDGMass();
  availableEnergy -= muonBindingEnergy;
  
  G4ThreeVector aMu3Mom(0,0,0);
  G4LorentzVector aMuMom(availableEnergy,aMu3Mom);
  
  // pick random proton.
  G4double residualMass=0;
  G4double targetMass=0;
  G4ThreeVector pResInCMS;
  G4double eEx=0;
  G4double eRest = 0;
  G4ReactionProduct* aNu=0;
  do
  {
    theN.Init(targetAtomicMass, targetCharge); 
    G4ThreeVector fermiMom;
    G4Nucleon * aNucleon = 0;
    G4int theProtonCounter = G4lrint( 0.5 + targetCharge * G4UniformRand() );
    G4int counter = 0;
    theN.StartLoop();
    while( (aNucleon=theN.GetNextNucleon()) )
    {
      if( aNucleon->GetDefinition() == G4Proton::ProtonDefinition() )
      {
	counter++;
	if(counter == theProtonCounter)
	{
          fermiMom = aNucleon->GetMomentum().vect();
	  break;
	}
      }
    }

    // Gett the nu momentum in the CMS
    G4double aNMass = G4Proton::ProtonDefinition()->GetPDGMass()/2.;
    G4LorentzVector theNeutronMom(std::sqrt(aNMass*aNMass+fermiMom.mag2()), fermiMom);
    G4LorentzVector theCMS = theNeutronMom+aMuMom;
    G4double p1 = (theCMS.mag()*theCMS.mag()-aNMass*aNMass)/(2.*theCMS.mag());
    G4LorentzRotation toCMS = theCMS.boostVector();
    G4LorentzRotation toLab = toCMS.inverse();

    // make the nu, and transform to lab;
    G4double cosTh = G4UniformRand();
    G4double phi = twopi*G4UniformRand();
    G4double theta = std::acos(cosTh);
    G4double sinth = std::sin(theta);
    G4ThreeVector randUnit(sinth*std::cos(phi), sinth*std::sin(phi), std::cos(theta) );
    G4LorentzVector finNuMom(p1, p1*randUnit);

    // make the neutrino an mu-neutrino with the above momentum and get the residual properties
    if(aNu) delete aNu;
    aNu = new G4ReactionProduct();
    aNu->SetDefinition( G4NeutrinoMu::NeutrinoMuDefinition() );
    pResInCMS = -finNuMom.vect();

    // nu in lab.
    finNuMom*=toLab;
    aNu->SetTotalEnergy( finNuMom.t() );
    aNu->SetMomentum( finNuMom.vect() );
    residualMass =  
             G4ParticleTable::GetParticleTable()
	     ->GetIonTable()->GetIonMass(G4lrint(targetCharge-1), G4lrint(targetAtomicMass));
    targetMass = 
             G4ParticleTable::GetParticleTable()
	     ->GetIonTable()->GetIonMass(G4lrint(targetCharge), G4lrint(targetAtomicMass));
    eRest = targetMass+availableEnergy-residualMass-finNuMom.t();

    // Call pre-compound on the rest.
    eEx = std::sqrt( (residualMass+eRest)*(residualMass+eRest) - pResInCMS.mag2() ) - residualMass;
  }
  while(eEx<=0);
  G4LorentzVector resV(residualMass+eRest, pResInCMS);
  G4LorentzRotation toBreit = resV.boostVector();
  G4LorentzRotation fromBreit = toBreit.inverse();
  resV*=toBreit;
  G4Fragment anInitialState;
  anInitialState.SetA(G4lrint(targetAtomicMass));
  anInitialState.SetZ(G4lrint(targetCharge-1));
  anInitialState.SetNumberOfParticles(2);
  anInitialState.SetNumberOfCharged(0);
  anInitialState.SetNumberOfHoles(1);
  anInitialState.SetMomentum(resV);
  G4ReactionProductVector * aPreResult = theHandler.BreakItUp(anInitialState);

  G4ReactionProductVector::iterator ires;
  G4double eBal = availableEnergy;
  for(ires=aPreResult->begin(); ires!=aPreResult->end(); ires++)
  {
    G4LorentzVector itV((*ires)->GetTotalEnergy(), (*ires)->GetMomentum());
    itV*=fromBreit;
    (*ires)->SetTotalEnergy(itV.t());
    (*ires)->SetMomentum(itV.vect());
    eBal-=itV.t()-itV.mag();
  }
  // fill neutrino into result
  aPreResult->push_back(aNu);
  eBal-=aNu->GetMomentum().mag();
  if(getenv("MuonMinusCaptureAtRestPrint"))
  {
    G4cout << "    Rough energy balance is "<<eBal<<" "
           <<availableEnergy<<" "<<aNu->GetMomentum().mag()<<G4endl;
  }    
  // return
  return aPreResult;
} 

