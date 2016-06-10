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
// $Id$
//
//   G4MuonMinusCaptureAtRest physics process
//   Larry Felawka (TRIUMF) and Art Olin (TRIUMF)
//   April 1998
//---------------------------------------------------------------------
//
// Modifications: 
// 18/08/2000  V.Ivanchenko Update description 
// 12/12/2003  H.P.Wellisch Completly rewrite mu-nuclear part
// 17/05/2006  V.Ivanchenko Cleanup
// 15/11/2006  V.Ivanchenko Review and rewrite all kinematics
// 24/01/2007  V.Ivanchenko Force to work with integer Z and A
// 23/01/2009  V.Ivanchenko Add deregistration
//
//-----------------------------------------------------------------------------

#include "G4MuonMinusCaptureAtRest.hh"
#include "G4DynamicParticle.hh"
#include "Randomize.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4He3.hh"
#include "G4NeutrinoMu.hh"
#include "G4Fragment.hh"
#include "G4ReactionProductVector.hh"
#include "G4Proton.hh"
#include "G4PionPlus.hh"
#include "G4GHEKinematicsVector.hh"
#include "G4Fancy3DNucleus.hh"
#include "G4ExcitationHandler.hh"
#include "G4HadronicProcessStore.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4MuonMinusCaptureAtRest::G4MuonMinusCaptureAtRest(const G4String& processName,
						   G4ProcessType   aType ) :
  G4VRestProcess (processName, aType), nCascade(0), targetZ(0), targetA(0), 
  isInitialised(false)
{
  SetProcessSubType(fHadronAtRest);
  Cascade    = new G4GHEKinematicsVector [17];
  pSelector  = new G4StopElementSelector();
  pEMCascade = new G4MuMinusCaptureCascade();
  theN       = new G4Fancy3DNucleus();
  theHandler = new G4ExcitationHandler();
  G4HadronicProcessStore::Instance()->RegisterExtraProcess(this);
  targetMass = 0.0;
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4MuonMinusCaptureAtRest::~G4MuonMinusCaptureAtRest()
{
  G4HadronicProcessStore::Instance()->DeRegisterExtraProcess(this);
  delete [] Cascade;
  delete pSelector;
  delete pEMCascade;
  delete theN;
  delete theHandler;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4MuonMinusCaptureAtRest::IsApplicable(const G4ParticleDefinition& p)
{
  return ( &p == G4MuonMinus::MuonMinus() );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4MuonMinusCaptureAtRest::PreparePhysicsTable(const G4ParticleDefinition& p) 
{
  G4HadronicProcessStore::Instance()->RegisterParticleForExtraProcess(this, &p);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4MuonMinusCaptureAtRest::BuildPhysicsTable(const G4ParticleDefinition& p) 
{
  G4HadronicProcessStore::Instance()->PrintInfo(&p);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* G4MuonMinusCaptureAtRest::AtRestDoIt(const G4Track& track, 
							const G4Step&)
{
  //
  // Handles MuonMinuss at rest; a MuonMinus can either create secondaries or
  // do nothing (in which case it should be sent back to decay-handling
  // section
  //
  aParticleChange.Initialize(track);

  // select element and get Z,A.
  G4Element* aEle = pSelector->GetElement(track.GetMaterial());
  targetZ = aEle->GetZ();
  targetA = G4double(G4int(aEle->GetN()+0.5)); 
  G4int ni = 0;

  G4IsotopeVector* isv = aEle->GetIsotopeVector();
  if(isv) ni = isv->size();

  if(ni == 1) {
    targetA = G4double(aEle->GetIsotope(0)->GetN());
  } else if(ni > 1) {
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
  nCascade   = 0;
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
  /*
  G4cout << "lambda= " << lambda << " lambdac= " << lambdac 
	 << " nem= " << nEmSecondaries << G4endl;
  */
  if( G4UniformRand()*lambda > lambdac) 
    pEMCascade->DoBoundMuonMinusDecay(targetZ, &nEmSecondaries, Cascade);
  else 
    captureResult = DoMuCapture();
  
  // fill the final state
  if(captureResult) nSecondaries += captureResult->size();
  else nSecondaries = nEmSecondaries;
  //G4cout << " nsec= " << nSecondaries << " nem= " << nEmSecondaries << G4endl;

  aParticleChange.SetNumberOfSecondaries( nSecondaries );

  G4double globalTime = track.GetGlobalTime();
  G4ThreeVector position = track.GetPosition();
  // Store nuclear cascade
  if(captureResult) {
    G4int n = captureResult->size();
    for ( G4int isec = 0; isec < n; isec++ ) {
      G4ReactionProduct* aParticle = captureResult->operator[](isec);
      G4DynamicParticle * aNewParticle = new G4DynamicParticle();
      aNewParticle->SetDefinition( aParticle->GetDefinition() );
      G4LorentzVector itV(aParticle->GetTotalEnergy(), aParticle->GetMomentum());
      aNewParticle->SetMomentum(itV.vect());
      G4double localtime = globalTime + tDelay + aParticle->GetTOF();
      G4Track* aNewTrack = new G4Track( aNewParticle, localtime, position);
      aNewTrack->SetTouchableHandle(track.GetTouchableHandle());
      aParticleChange.AddSecondary( aNewTrack );
      delete aParticle;
    }
    delete captureResult;
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ReactionProductVector* G4MuonMinusCaptureAtRest::DoMuCapture()
{
  G4double mumass = G4MuonMinus::MuonMinus()->GetPDGMass();
  G4double muBindingEnergy = pEMCascade->GetKShellEnergy(targetZ);
  /*
  G4cout << "G4MuonMinusCaptureAtRest::DoMuCapture called Emu= "
	 << muBindingEnergy << G4endl;
  */
  // Energy on K-shell
  G4double muEnergy = mumass + muBindingEnergy;
  G4double muMom = std::sqrt(muBindingEnergy*(muBindingEnergy + 2.0*mumass));
  G4double availableEnergy = targetMass + mumass - muBindingEnergy;
  G4LorentzVector momInitial(0.0,0.0,0.0,availableEnergy);
  G4LorentzVector momResidual;

  G4ThreeVector vmu = muMom*pEMCascade->GetRandomVec();
  G4LorentzVector aMuMom(vmu, muEnergy);

  G4double residualMass = 
    G4NucleiProperties::GetNuclearMass(targetA, targetZ - 1.0);

  G4ReactionProductVector* aPreResult = 0;
  G4ReactionProduct* aNu = new G4ReactionProduct();
  aNu->SetDefinition( G4NeutrinoMu::NeutrinoMu() );

  G4int iz = G4int(targetZ);
  G4int ia = G4int(targetA);

  // p, d, t, 3He or alpha as target
  if(iz <= 2) {

    if(ia > 1) {
      if(iz == 1 && ia == 2) { 
	availableEnergy -= neutron_mass_c2;
      } else if(iz == 1 && ia == 3) {
	availableEnergy -= 2.0*neutron_mass_c2;
      } else if(iz == 2) {
        G4ParticleDefinition* pd = 0;
	if (ia == 3) {
          pd = G4Deuteron::Deuteron();
        } else if(ia == 4) {
          pd = G4Triton::Triton();
        } else { 
	  pd = G4ParticleTable::GetParticleTable()->FindIon(1,ia-1,0,1);
        }

	//	G4cout << "Extra " << pd->GetParticleName() << G4endl;
	availableEnergy -= pd->GetPDGMass();
      }
    }
    //
    //  Computation in assumption of CM collision of mu and nucleaon
    //  
    G4double Enu  = 0.5*(availableEnergy - 
			 neutron_mass_c2*neutron_mass_c2/availableEnergy);

    // make the nu, and transform to lab;
    G4ThreeVector nu3Mom = Enu*pEMCascade->GetRandomVec();

    G4ReactionProduct* aN = new G4ReactionProduct();
    aN->SetDefinition( G4Neutron::Neutron() );
    aN->SetTotalEnergy( availableEnergy - Enu );
    aN->SetMomentum( -nu3Mom );

    aNu->SetTotalEnergy( Enu );
    aNu->SetMomentum( nu3Mom );
    aPreResult = new G4ReactionProductVector();

    aPreResult->push_back(aN ); 
    aPreResult->push_back(aNu);

    if(verboseLevel > 1)
      G4cout << "DoMuCapture on H or He" 
	     <<" EkinN(MeV)= " << (availableEnergy - Enu - neutron_mass_c2)/MeV
	     <<" Enu(MeV)= "<<aNu->GetTotalEnergy()/MeV
	     <<" n= " << aPreResult->size()
	     <<G4endl;

    return aPreResult;
  }

  // pick random proton inside nucleus 
  G4double eEx;
  do {
    theN->Init(ia, iz); 
    G4LorentzVector thePMom;
    G4Nucleon * aNucleon = 0;
    G4int theProtonCounter = G4int( targetZ * G4UniformRand() );
    G4int counter = 0;
    theN->StartLoop();

    while( (aNucleon=theN->GetNextNucleon()) ) {

      if( aNucleon->GetDefinition() == G4Proton::Proton() ) {
	counter++;
	if(counter == theProtonCounter) {
	  thePMom  = aNucleon->GetMomentum();
	  break;
	}
      }
    }

    // Get the nu momentum in the CMS
    G4LorentzVector theCMS = thePMom + aMuMom;
    G4ThreeVector bst = theCMS.boostVector();

    G4double Ecms = theCMS.mag();
    G4double Enu  = 0.5*(Ecms - neutron_mass_c2*neutron_mass_c2/Ecms);
    eEx = 0.0;

    if(Enu > 0.0) {
      // make the nu, and transform to lab;
      G4ThreeVector nu3Mom = Enu*pEMCascade->GetRandomVec();
      G4LorentzVector nuMom(nu3Mom, Enu);

      // nu in lab.
      nuMom.boost(bst);
      aNu->SetTotalEnergy( nuMom.e() );
      aNu->SetMomentum( nuMom.vect() );
    
      // make residual
      momResidual = momInitial - nuMom;

      // Call pre-compound on the rest.
      eEx = momResidual.mag();
      if(verboseLevel > 1)
	G4cout << "G4MuonMinusCaptureAtRest::DoMuCapture: " 
	       << " Eex(MeV)= " << (eEx-residualMass)/MeV
	       << " Enu(MeV)= "<<aNu->GetTotalEnergy()/MeV
	       <<G4endl;
    }
  } while(eEx <= residualMass);

//  G4cout << "muonCapture : " << eEx << " " << residualMass 
//         << " A,Z= " << targetA << ", "<< targetZ 
//	 << "  " << G4int(targetA) << ", " << G4int(targetZ) << G4endl;

  //
  // Start Deexcitation
  //
  G4ThreeVector fromBreit = momResidual.boostVector();
  G4LorentzVector fscm(0.0,0.0,0.0, eEx);
  G4Fragment anInitialState;
  anInitialState.SetA(targetA);
  anInitialState.SetZ(G4double(iz - 1));
  anInitialState.SetNumberOfParticles(2);
  anInitialState.SetNumberOfCharged(0);
  anInitialState.SetNumberOfHoles(1);
  anInitialState.SetMomentum(fscm);
  aPreResult = theHandler->BreakItUp(anInitialState);

  G4ReactionProductVector::iterator ires;
  G4double eBal = availableEnergy - aNu->GetTotalEnergy();
  for(ires=aPreResult->begin(); ires!=aPreResult->end(); ires++) {
    G4LorentzVector itV((*ires)->GetTotalEnergy(), (*ires)->GetMomentum());
    itV.boost(fromBreit);
    (*ires)->SetTotalEnergy(itV.t());
    (*ires)->SetMomentum(itV.vect());
    eBal -= itV.t();
  }
  //
  // fill neutrino into result
  //
  aPreResult->push_back(aNu);
 
  if(verboseLevel > 1)
    G4cout << "DoMuCapture:  Nsec= " 
	   << aPreResult->size() << " Ebalance(MeV)= " << eBal/MeV
	   <<" E0(MeV)= " <<availableEnergy/MeV
	   <<" Mres(GeV)= " <<residualMass/GeV
	   <<G4endl;

  return aPreResult;
} 

