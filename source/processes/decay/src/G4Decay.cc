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
// $Id: G4Decay.cc,v 1.15 2002-02-12 08:57:20 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
#include "G4DynamicParticle.hh"
#include "G4DecayProducts.hh"
#include "G4DecayTable.hh"
#include "G4PhysicsLogVector.hh"
#include "G4ParticleChangeForDecay.hh"
#include "G4VExtDecayer.hh"

// constructor
G4Decay::G4Decay(const G4String& processName)
                               :G4VRestDiscreteProcess(processName, fDecay),
				verboseLevel(1),
                                HighestValue(20.0),
                                pExtDecayer(0)
{
#ifdef G4VERBOSE
  if (GetVerboseLevel()>1) {
    G4cerr << "G4Decay  constructor " << "  Name:" << processName << G4endl;
  }
#endif
  pParticleChange = &fParticleChangeForDecay;
}

G4Decay::~G4Decay()
{
  if (pExtDecayer) {
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

G4double G4Decay::GetMeanLifeTime(const G4Track&    aTrack,
                                  G4ForceCondition*)
{
   // get particle type 
   const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();

   // returns the mean free path in GEANT4 internal units
   G4double meanlife;
   G4ParticleDefinition* aParticleDef = aParticle->GetDefinition();
   G4double aLife = aParticleDef->GetPDGLifeTime();

   // check if the particle is stable?
   if (aParticleDef->GetPDGStable()) {
     meanlife = DBL_MAX;
    
   } else if (aLife < 0.0) {
     meanlife = DBL_MAX;
      
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

   // returns the mean free path in GEANT4 internal units
   G4double pathlength;
   G4ParticleDefinition* aParticleDef = aParticle->GetDefinition();
   G4double aCtau = c_light * aParticleDef->GetPDGLifeTime();
   G4double aMass = aParticle->GetMass();

   // check if the particle is stable?
   if (aParticleDef->GetPDGStable()) {
     pathlength = DBL_MAX;

   } else if (aCtau < 0.0) {
     pathlength =  DBL_MAX;
      
   //check if the particle has very short life time ?
   } else if (aCtau < DBL_MIN) { 
     pathlength =  DBL_MIN;
 
   //check if zero mass
   } else if (aMass <  DBL_MIN)  {
     pathlength =  DBL_MAX;

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
	 G4cerr << "G4Decay::GetMeanFreePath()   !!particle stops!!";
         G4cerr << aParticleDef->GetParticleName() << G4endl;
	 G4cerr << "KineticEnergy:" << aParticle->GetKineticEnergy()/GeV <<"[GeV]";
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
  G4ParticleDefinition* aParticleDef = aParticle->GetDefinition();

  // check if  the particle is stable
  if (aParticleDef->GetPDGStable()) return &fParticleChangeForDecay ;
 

  //check if thePreAssignedDecayProducts exists
  const G4DecayProducts* o_products = (aParticle->GetPreAssignedDecayProducts());
  G4bool isPreAssigned = (o_products != 0);   
  G4DecayProducts* products = 0;

  // decay table
  G4DecayTable   *decaytable = aParticleDef->GetDecayTable();
 
  // check if external decayer exists
  G4bool isExtDecayer = (decaytable == 0) && (pExtDecayer !=0);

  // Error due to NO Decay Table 
  if ( (decaytable == 0) && !isExtDecayer &&!isPreAssigned ){
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cerr <<  "G4Decay::DoIt  : decay table not defined  for ";
      G4cerr << aParticle->GetDefinition()->GetParticleName()<< G4endl;
    }
#endif
    fParticleChangeForDecay.SetNumberOfSecondaries(0);
    // Kill the parent particle
    fParticleChangeForDecay.SetStatusChange( fStopAndKill ) ;
    fParticleChangeForDecay.SetLocalEnergyDeposit(0.0); 
    
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
    // decay acoording to decay table
    // choose a decay channel
    G4VDecayChannel *decaychannel = decaytable->SelectADecayChannel();
    if (decaychannel == 0 ){
      // decay channel not found
      G4Exception("G4Decay::DoIt  : can not determine decay channel ");
    } else {
      G4int temp = decaychannel->GetVerboseLevel();
      // execute DecayIt() 
#ifdef G4VERBOSE
      if (GetVerboseLevel()>1) {
	G4cerr << "G4Decay::DoIt  : selected decay channel  addr:" << decaychannel <<G4endl;
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
      // for debug
      //if (! products->IsChecked() ) products->DumpInfo();
#endif
    }
  }
  
  // get parent particle information ...................................
  G4double   ParentEnergy  = aParticle->GetTotalEnergy();
  G4ThreeVector ParentDirection(aParticle->GetMomentumDirection());

  //boost all decay products to laboratory frame
  G4double energyDeposit = 0.0;
  G4double finalGlobalTime = aTrack.GetGlobalTime();
  if (aTrack.GetTrackStatus() == fStopButAlive ){
    // AtRest case
    finalGlobalTime += fRemainderLifeTime;
    energyDeposit += aParticle->GetKineticEnergy();
    if (isPreAssigned) products->Boost( ParentEnergy, ParentDirection); 
  } else {
    // PostStep case
    if (!isExtDecayer) products->Boost( ParentEnergy, ParentDirection); 
  }

  //add products in fParticleChangeForDecay
  G4int numberOfSecondaries = products->entries();
  fParticleChangeForDecay.SetNumberOfSecondaries(numberOfSecondaries);
#ifdef G4VERBOSE
  if (GetVerboseLevel()>1) {
    G4cerr << "G4Decay::DoIt  : Decay vertex :";
    G4cerr << " Time: " << finalGlobalTime/ns << "[ns]";
    G4cerr << " X:" << (aTrack.GetPosition()).x() /cm << "[cm]";
    G4cerr << " Y:" << (aTrack.GetPosition()).y() /cm << "[cm]";
    G4cerr << " Z:" << (aTrack.GetPosition()).z() /cm << "[cm]";
    G4cerr << G4endl;
    G4cerr << "G4Decay::DoIt  : decay products in Lab. Frame" << G4endl;
    products->DumpInfo();
  }
#endif
  G4int index;
  G4ThreeVector currentPosition;
  for (index=0; index < numberOfSecondaries; index++) 
  { 
     // get current position of the track
     currentPosition = aTrack.GetPosition();
     // create a new track object
     G4Track* secondary = new G4Track( products->PopProducts(),
				      finalGlobalTime ,
				      currentPosition );
     // switch on good for tracking flag 
     secondary->SetGoodForTrackingFlag();
     // add the secondary track in the List
     fParticleChangeForDecay.AddSecondary(secondary);
  }
  delete products;

  // Kill the parent particle
  fParticleChangeForDecay.SetStatusChange( fStopAndKill ) ;
  fParticleChangeForDecay.SetLocalEnergyDeposit(energyDeposit); 
  fParticleChangeForDecay.SetTimeChange( finalGlobalTime );
  // reset NumberOfInteractionLengthLeft
  ClearNumberOfInteractionLengthLeft();

  return &fParticleChangeForDecay ;
} 






