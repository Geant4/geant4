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
// $Id: G4Decay.cc 105727 2017-08-16 12:47:05Z gcosmo $
//
// 
// --------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      History: first implementation, based on object model of
//      2nd December 1995, G.Cosmo
//      7 July 1996 H.Kurashige
// ------------------------------------------------------------
//   remove  BuildPhysicsTable()  28 Nov. 1997  H.Kurashige
//   change  DBL_EPSIRON to DBL_MIN 14 Dec. 1997  H.Kurashige
//   modified for new ParticleChange 12 Mar. 1998  H.Kurashige
//   modified for "GoodForTrackingFlag" 19 June 1998  H.Kurashige
//   rename thePhysicsTable to aPhyscisTable 2 Aug. 1998 H.Kurashige
//   modified IsApplicable in order to protect the decay from registered 
//   to resonances    12 Dec. 1998   H.Kurashige 
//   remove G4ParticleMomentum  6 Feb. 99 H.Kurashige
//   modified  IsApplicable to activate G4Decay for resonances  1 Mar. 00 H.Kurashige 
//   Add External Decayer         23 Feb. 2001  H.Kurashige
//   change LowestBinValue,HighestBinValue and TotBin(200) 9 Feb. 2002
//

#include "G4Decay.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4DynamicParticle.hh"
#include "G4DecayProducts.hh"
#include "G4DecayTable.hh"
#include "G4VDecayChannel.hh"
#include "G4PhysicsLogVector.hh"
#include "G4ParticleChangeForDecay.hh"
#include "G4VExtDecayer.hh"

// constructor
G4Decay::G4Decay(const G4String& processName)
                               :G4VRestDiscreteProcess(processName, fDecay),
				verboseLevel(1),
                                HighestValue(20.0),
				fRemainderLifeTime(-1.0),
                                pExtDecayer(nullptr)
{
  // set Process Sub Type
  SetProcessSubType(static_cast<int>(DECAY));

#ifdef G4VERBOSE
  if (GetVerboseLevel()>1) {
    G4cout << "G4Decay  constructor " << "  Name:" << processName << G4endl;
  }
#endif

  pParticleChange = &fParticleChangeForDecay;
}

G4Decay::~G4Decay()
{
  if (pExtDecayer != nullptr) {
    delete pExtDecayer;
  }
}

G4bool G4Decay::IsApplicable(const G4ParticleDefinition& aParticleType)
{
   // check if the particle is stable?
   if (aParticleType.GetPDGLifeTime() <0.0) {
     return false;
   } else if (aParticleType.GetPDGMass() <= 0.0*MeV) {
     return false;
   } else {
     return true; 
   }
}

G4double G4Decay::GetMeanLifeTime(const G4Track& aTrack  ,
                                  G4ForceCondition*)
{
   // returns the mean free path in GEANT4 internal units
   G4double meanlife;

   // get particle 
   const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
   const G4ParticleDefinition* aParticleDef = aParticle->GetDefinition();
   G4double aLife = aParticleDef->GetPDGLifeTime();

   // check if the particle is stable?
   if (aParticleDef->GetPDGStable()) {
    //1000000 times the life time of the universe
     meanlife = 1e24 * s;
    
   } else {
     meanlife = aLife;
   }

#ifdef G4VERBOSE
   if (GetVerboseLevel()>1) {
     G4cout << "mean life time: "<< meanlife/ns << "[ns]" << G4endl;
   }
#endif

   return  meanlife;
}

G4double G4Decay::GetMeanFreePath(const G4Track& aTrack,G4double, G4ForceCondition*)
{
   // get particle 
   const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
   const G4ParticleDefinition* aParticleDef = aParticle->GetDefinition();
   G4double aMass = aParticle->GetMass();
   G4double aLife = aParticleDef->GetPDGLifeTime();


    // returns the mean free path in GEANT4 internal units
   G4double pathlength;
   G4double aCtau = c_light * aLife;

   // check if the particle is stable?
   if (aParticleDef->GetPDGStable()) {
     pathlength = DBL_MAX;

   //check if the particle has very short life time ?
   } else if (aCtau < DBL_MIN) { 
     pathlength =  DBL_MIN;
 
   } else {
    //calculate the mean free path 
    // by using normalized kinetic energy (= Ekin/mass)
     G4double   rKineticEnergy = aParticle->GetKineticEnergy()/aMass; 
     if ( rKineticEnergy > HighestValue) {
       // gamma >>  1
       pathlength = ( rKineticEnergy + 1.0)* aCtau;
     } else if ( rKineticEnergy < DBL_MIN ) {
       // too slow particle
#ifdef G4VERBOSE
       if (GetVerboseLevel()>1) {
	 G4cout << "G4Decay::GetMeanFreePath()   !!particle stops!!";
         G4cout << aParticleDef->GetParticleName() << G4endl;
	 G4cout << "KineticEnergy:" << aParticle->GetKineticEnergy()/GeV <<"[GeV]";
       }
#endif
       pathlength = DBL_MIN;
     } else {
       // beta <1 
       pathlength = (aParticle->GetTotalMomentum())/aMass*aCtau ;
     }
   }
  return  pathlength;
}

void G4Decay::BuildPhysicsTable(const G4ParticleDefinition&)
{
  return;
}

G4VParticleChange* G4Decay::DecayIt(const G4Track& aTrack, const G4Step& )
{
  // The DecayIt() method returns by pointer a particle-change object.
  // Units are expressed in GEANT4 internal units.

  //   Initialize ParticleChange
  //     all members of G4VParticleChange are set to equal to 
  //     corresponding member in G4Track
  fParticleChangeForDecay.Initialize(aTrack);

  // get particle 
  const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
  const G4ParticleDefinition* aParticleDef = aParticle->GetDefinition();

  // check if  the particle is stable
  if (aParticleDef->GetPDGStable()) return &fParticleChangeForDecay ;
 

  //check if thePreAssignedDecayProducts exists
  const G4DecayProducts* o_products = (aParticle->GetPreAssignedDecayProducts());
  G4bool isPreAssigned = (o_products != nullptr);   
  G4DecayProducts* products = nullptr;

  // decay table
  G4DecayTable   *decaytable = aParticleDef->GetDecayTable();
 
  // check if external decayer exists
  G4bool isExtDecayer = (decaytable == nullptr) && (pExtDecayer != nullptr);

  // Error due to NO Decay Table 
  if ( (decaytable == nullptr) && !isExtDecayer && !isPreAssigned ){
    if (GetVerboseLevel()>0) {
      G4cout <<  "G4Decay::DoIt  : decay table not defined  for ";
      G4cout << aParticle->GetDefinition()->GetParticleName()<< G4endl;
    }
    G4Exception( "G4Decay::DecayIt ",
                 "DECAY101",JustWarning, 
                 "Decay table is not defined");

    fParticleChangeForDecay.SetNumberOfSecondaries(0);
    // Kill the parent particle
    fParticleChangeForDecay.ProposeTrackStatus( fStopAndKill ) ;
    fParticleChangeForDecay.ProposeLocalEnergyDeposit(0.0); 
    
    ClearNumberOfInteractionLengthLeft();
    return &fParticleChangeForDecay ;
  }

  if (isPreAssigned) {
    // copy decay products 
    products = new G4DecayProducts(*o_products); 
  } else if ( isExtDecayer ) {
    // decay according to external decayer
    products = pExtDecayer->ImportDecayProducts(aTrack);
  } else {
    // Decay according to decay table.
    // Keep trying to choose a candidate decay channel if the dynamic mass
    // of the decaying particle is below the sum of the PDG masses of the
    // candidate daughter particles.
    // This is needed because the decay table used in Geant4 is based on
    // the assumption of nominal PDG masses, but a wide resonance can have
    // a dynamic masses well below its nominal PDG masses, and therefore
    // some of its decay channels can be below the kinematical threshold. 
    // Note that, for simplicity, we ignore here the possibility that
    // one or more of the candidate daughter particles can be, in turn,
    // wide resonance. However, if this is the case, and the channel is
    // accepted, then the masses of the resonance daughter particles will
    // be sampled by taking into account their widths.
    G4VDecayChannel* decaychannel = nullptr;
    G4double massParent = aParticle->GetMass();
    decaychannel = decaytable->SelectADecayChannel(massParent);
    if ( decaychannel == nullptr) {
      // decay channel not found
	   G4ExceptionDescription ed;
      ed << "Can not determine decay channel for " 
	 << aParticleDef->GetParticleName() << G4endl 
	 << "  mass of dynamic particle: " 
	 << massParent/GeV << " (GEV)" << G4endl
	 << "  dacay table has " << decaytable->entries() 
	 << " entries" << G4endl;
      G4double checkedmass=massParent;
      if (massParent < 0.) {
	checkedmass=aParticleDef->GetPDGMass();
	ed << "Using PDG mass ("<<checkedmass/GeV 
	   << "(GeV)) in IsOKWithParentMass" << G4endl;	
      }
      for (G4int ic =0;ic <decaytable->entries();++ic) {
	G4VDecayChannel * dc= decaytable->GetDecayChannel(ic);
	ed << ic << ": BR " << dc->GetBR() << ", IsOK? " 
	   << dc->IsOKWithParentMass(checkedmass)
	   << ", --> "; 
	G4int ndaughters=dc->GetNumberOfDaughters();
	for (G4int id=0;id<ndaughters;++id) {
	  if (id>0) ed << " + ";   // seperator, except for first
	  ed << dc->GetDaughterName(id);
	}
	ed << G4endl;
      }
      G4Exception("G4Decay::DoIt", "DECAY003", FatalException,ed);
    } else {
      // execute DecayIt() 
#ifdef G4VERBOSE
      G4int temp = decaychannel->GetVerboseLevel();
      if (GetVerboseLevel()>1) {
	G4cout << "G4Decay::DoIt  : selected decay channel  addr:" 
	       << decaychannel <<G4endl;
	decaychannel->SetVerboseLevel(GetVerboseLevel());
      }
#endif
      products = decaychannel->DecayIt(aParticle->GetMass());
#ifdef G4VERBOSE
      if (GetVerboseLevel()>1) {
	decaychannel->SetVerboseLevel(temp);
      }
#endif
#ifdef G4VERBOSE
      if (GetVerboseLevel()>2) {
	if (! products->IsChecked() ) products->DumpInfo();
      }
#endif
    }
  }
  
  // get parent particle information ...................................
  G4double   ParentEnergy  = aParticle->GetTotalEnergy();
  G4double   ParentMass    = aParticle->GetMass();
  if (ParentEnergy < ParentMass) {
    if (GetVerboseLevel()>0) {
      G4cout << "G4Decay::DoIt  : Total Energy is less than its mass" << G4endl;
      G4cout << " Particle: " << aParticle->GetDefinition()->GetParticleName();
      G4cout << " Energy:"    << ParentEnergy/MeV << "[MeV]";
      G4cout << " Mass:"    << ParentMass/MeV << "[MeV]";
      G4cout << G4endl;
    }
    G4Exception( "G4Decay::DecayIt ",
                 "DECAY102",JustWarning, 
                 "Total Energy is less than its mass");
    ParentEnergy = ParentMass;
  }

  G4ThreeVector ParentDirection(aParticle->GetMomentumDirection());

  //boost all decay products to laboratory frame
  G4double energyDeposit = 0.0;
  G4double finalGlobalTime = aTrack.GetGlobalTime();
  G4double finalLocalTime = aTrack.GetLocalTime();
  if (aTrack.GetTrackStatus() == fStopButAlive ){
    // AtRest case
    finalGlobalTime += fRemainderLifeTime;
    finalLocalTime += fRemainderLifeTime;
    energyDeposit += aParticle->GetKineticEnergy();
    if (isPreAssigned) products->Boost( ParentEnergy, ParentDirection);
  } else {
    // PostStep case
    if (!isExtDecayer) products->Boost( ParentEnergy, ParentDirection);
  }
   // set polarization for daughter particles
   DaughterPolarization(aTrack, products);


  //add products in fParticleChangeForDecay
  G4int numberOfSecondaries = products->entries();
  fParticleChangeForDecay.SetNumberOfSecondaries(numberOfSecondaries);
#ifdef G4VERBOSE
  if (GetVerboseLevel()>1) {
    G4cout << "G4Decay::DoIt  : Decay vertex :";
    G4cout << " Time: " << finalGlobalTime/ns << "[ns]";
    G4cout << " X:" << (aTrack.GetPosition()).x() /cm << "[cm]";
    G4cout << " Y:" << (aTrack.GetPosition()).y() /cm << "[cm]";
    G4cout << " Z:" << (aTrack.GetPosition()).z() /cm << "[cm]";
    G4cout << G4endl;
    G4cout << "G4Decay::DoIt  : decay products in Lab. Frame" << G4endl;
    products->DumpInfo();
  }
#endif
  G4int index;
  G4ThreeVector currentPosition;
  const G4TouchableHandle thand = aTrack.GetTouchableHandle();
  for (index=0; index < numberOfSecondaries; index++){
    // get current position of the track
    currentPosition = aTrack.GetPosition();
    // create a new track object
    G4Track* secondary = new G4Track( products->PopProducts(),
				      finalGlobalTime ,
				      currentPosition );
    // switch on good for tracking flag
    secondary->SetGoodForTrackingFlag();
    secondary->SetTouchableHandle(thand);
    // add the secondary track in the List
    fParticleChangeForDecay.AddSecondary(secondary);
  }
  delete products;
  
  // Kill the parent particle
  fParticleChangeForDecay.ProposeTrackStatus( fStopAndKill ) ;
  fParticleChangeForDecay.ProposeLocalEnergyDeposit(energyDeposit); 
  fParticleChangeForDecay.ProposeLocalTime( finalLocalTime );
  
  // Clear NumberOfInteractionLengthLeft
  ClearNumberOfInteractionLengthLeft();
  
  return &fParticleChangeForDecay ;
} 

void G4Decay::DaughterPolarization(const G4Track& , G4DecayProducts* )
{
  // empty implementation
}



void G4Decay::StartTracking(G4Track*)
{
  currentInteractionLength = -1.0;
  ResetNumberOfInteractionLengthLeft();
 
  fRemainderLifeTime = -1.0;
}

void G4Decay::EndTracking()
{
  // Clear NumberOfInteractionLengthLeft
  ClearNumberOfInteractionLengthLeft();

  currentInteractionLength = -1.0;
}


G4double G4Decay::PostStepGetPhysicalInteractionLength(
                             const G4Track& track,
                             G4double   previousStepSize,
                             G4ForceCondition* condition
                            )
{
   // condition is set to "Not Forced"
  *condition = NotForced;

   // pre-assigned Decay time
  G4double pTime = track.GetDynamicParticle()->GetPreAssignedDecayProperTime();
  G4double aLife = track.GetDynamicParticle()->GetDefinition()->GetPDGLifeTime();

  if (pTime < 0.) {
    // normal case 
    if ( previousStepSize > 0.0){
      // subtract NumberOfInteractionLengthLeft 
      SubtractNumberOfInteractionLengthLeft(previousStepSize);
      if(theNumberOfInteractionLengthLeft<0.){
	theNumberOfInteractionLengthLeft=perMillion;
      }
      fRemainderLifeTime = theNumberOfInteractionLengthLeft*aLife;
    }
    // get mean free path
    currentInteractionLength = GetMeanFreePath(track, previousStepSize, condition);
    
#ifdef G4VERBOSE
    if ((currentInteractionLength <=0.0) || (verboseLevel>2)){
      G4cout << "G4Decay::PostStepGetPhysicalInteractionLength " << G4endl;
      track.GetDynamicParticle()->DumpInfo();
      G4cout << " in Material  " << track.GetMaterial()->GetName() <<G4endl;
      G4cout << "MeanFreePath = " << currentInteractionLength/cm << "[cm]" <<G4endl;
    }
#endif

    G4double value;
    if (currentInteractionLength <DBL_MAX) {
      value = theNumberOfInteractionLengthLeft * currentInteractionLength;
      //fRemainderLifeTime = theNumberOfInteractionLengthLeft*aLife;
    } else {
      value = DBL_MAX;
    }

    return value;

  } else {
    //pre-assigned Decay time case
    // reminder proper time
    fRemainderLifeTime = pTime - track.GetProperTime();
    if (fRemainderLifeTime <= 0.0) fRemainderLifeTime = DBL_MIN;
    
    G4double  rvalue=0.0; 
    // use pre-assigned Decay time to determine PIL
    if (aLife>0.0) {
      // ordinary particle
      rvalue = (fRemainderLifeTime/aLife)*GetMeanFreePath(track, previousStepSize, condition);
    } else {
     // shortlived particle
      rvalue = c_light * fRemainderLifeTime;
     // by using normalized kinetic energy (= Ekin/mass)
     G4double   aMass =  track.GetDynamicParticle()->GetMass();
     rvalue *= track.GetDynamicParticle()->GetTotalMomentum()/aMass;
    }
    return rvalue;
  }
}

G4double G4Decay::AtRestGetPhysicalInteractionLength(
                             const G4Track& track,
                             G4ForceCondition* condition
                            )
{
     // condition is set to "Not Forced"
  *condition = NotForced;

  G4double pTime = track.GetDynamicParticle()->GetPreAssignedDecayProperTime();
  if (pTime >= 0.) {
    fRemainderLifeTime = pTime - track.GetProperTime();
    if (fRemainderLifeTime <= 0.0) fRemainderLifeTime = DBL_MIN;
  } else {
    fRemainderLifeTime = 
      theNumberOfInteractionLengthLeft * GetMeanLifeTime(track, condition);
  }
  return fRemainderLifeTime;
}


void G4Decay::SetExtDecayer(G4VExtDecayer* val)
{
  pExtDecayer = val;

  // set Process Sub Type
  if ( pExtDecayer !=0 ) {
    SetProcessSubType(static_cast<int>(DECAY_External));
  }
}

G4VParticleChange* G4Decay::PostStepDoIt(
			     const G4Track& aTrack,
			     const G4Step&  aStep
			    )
{
  if ( (aTrack.GetTrackStatus() == fStopButAlive ) ||
       (aTrack.GetTrackStatus() == fStopAndKill )   ){
    fParticleChangeForDecay.Initialize(aTrack);
    return &fParticleChangeForDecay;
  } else {
    return DecayIt(aTrack, aStep);
  }
}

void G4Decay::ProcessDescription(std::ostream& outFile) const
{
  outFile << GetProcessName() << ": Decay of particles. \n"
	  << "kinematics of daughters are dertermined by DecayChannels " 
          << " or by PreAssignedDecayProducts\n";
}
