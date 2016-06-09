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
// $Id: G4UnknownDecay.cc,v 1.2 2004/10/19 00:56:34 kurasige Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// 
// --------------------------------------------------------------
//      GEANT 4 class implementation file
//
// ------------------------------------------------------------
//

#include "G4UnknownDecay.hh"
#include "G4DynamicParticle.hh"
#include "G4DecayProducts.hh"
#include "G4PhysicsLogVector.hh"
#include "G4ParticleChangeForDecay.hh"

// constructor
G4UnknownDecay::G4UnknownDecay(const G4String& processName)
                               :G4VDiscreteProcess(processName, fDecay),
				verboseLevel(1),
                                HighestValue(20.0)
{
#ifdef G4VERBOSE
  if (GetVerboseLevel()>1) {
    G4cerr << "G4UnknownDecay  constructor " << "  Name:" << processName << G4endl;
  }
#endif
  pParticleChange = &fParticleChangeForDecay;
}

G4UnknownDecay::~G4UnknownDecay()
{
}

G4bool G4UnknownDecay::IsApplicable(const G4ParticleDefinition& aParticleType)
{
  if(aParticleType.GetParticleName()=="unknown") return true;
  return false;
}

G4double G4UnknownDecay::GetMeanFreePath(const G4Track& /*aTrack*/,G4double, G4ForceCondition*)
{
//   // get particle 
//   const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
//
//   // returns the mean free path in GEANT4 internal units
//   G4double pathlength;
//   G4double pTime = aParticle->GetPreAssignedDecayProperTime();
//   if(pTime<0.) pTime = DBL_MIN;
//   G4double aCtau = c_light * pTime;
//   G4double aMass = aParticle->GetMass();
//
//   G4double rKineticEnergy = aParticle->GetKineticEnergy()/aMass; 
//   if ( rKineticEnergy > HighestValue) {
//     // gamma >>  1
//     pathlength = ( rKineticEnergy + 1.0)* aCtau;
//   } else if ( rKineticEnergy < DBL_MIN ) {
//     // too slow particle
//     pathlength = DBL_MIN;
//   } else {
//     pathlength = (aParticle->GetTotalMomentum())/aMass*aCtau ;
//   }
//   return  pathlength;
   return  DBL_MIN;
}

void G4UnknownDecay::BuildPhysicsTable(const G4ParticleDefinition&)
{
  return;
}

G4VParticleChange* G4UnknownDecay::DecayIt(const G4Track& aTrack, const G4Step& )
{
  // The DecayIt() method returns by pointer a particle-change object.
  // Units are expressed in GEANT4 internal units.

  //   Initialize ParticleChange
  //     all members of G4VParticleChange are set to equal to 
  //     corresponding member in G4Track
  fParticleChangeForDecay.Initialize(aTrack);

  // get particle 
  const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();

  //check if thePreAssignedDecayProducts exists
  const G4DecayProducts* o_products = (aParticle->GetPreAssignedDecayProducts());
  G4bool isPreAssigned = (o_products != 0);   
  G4DecayProducts* products = 0;

  if (!isPreAssigned ){
    fParticleChangeForDecay.SetNumberOfSecondaries(0);
    // Kill the parent particle
    fParticleChangeForDecay.ProposeTrackStatus( fStopAndKill ) ;
    fParticleChangeForDecay.ProposeLocalEnergyDeposit(0.0); 
    
    ClearNumberOfInteractionLengthLeft();
    return &fParticleChangeForDecay ;
  }

  // copy decay products 
  products = new G4DecayProducts(*o_products); 
  
  // get parent particle information ...................................
  G4double   ParentEnergy  = aParticle->GetTotalEnergy();
  G4ThreeVector ParentDirection(aParticle->GetMomentumDirection());

  G4double energyDeposit = 0.0;
  G4double finalGlobalTime = aTrack.GetGlobalTime();
  //boost all decay products to laboratory frame
  //if the particle has traveled 
  if(aParticle->GetPreAssignedDecayProperTime()>0.) {
    products->Boost( ParentEnergy, ParentDirection);
  }

  //add products in fParticleChangeForDecay
  G4int numberOfSecondaries = products->entries();
  fParticleChangeForDecay.SetNumberOfSecondaries(numberOfSecondaries);
#ifdef G4VERBOSE
  if (GetVerboseLevel()>1) {
    G4cerr << "G4UnknownDecay::DoIt  : Decay vertex :";
    G4cerr << " Time: " << finalGlobalTime/ns << "[ns]";
    G4cerr << " X:" << (aTrack.GetPosition()).x() /cm << "[cm]";
    G4cerr << " Y:" << (aTrack.GetPosition()).y() /cm << "[cm]";
    G4cerr << " Z:" << (aTrack.GetPosition()).z() /cm << "[cm]";
    G4cerr << G4endl;
    G4cerr << "G4UnknownDecay::DoIt  : decay products in Lab. Frame" << G4endl;
    products->DumpInfo();
  }
#endif
  G4int index;
  G4ThreeVector currentPosition;
  const G4TouchableHandle thand = aTrack.GetTouchableHandle();
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
     secondary->SetTouchableHandle(thand);
     // add the secondary track in the List
     fParticleChangeForDecay.AddSecondary(secondary);
  }
  delete products;

  // Kill the parent particle
  fParticleChangeForDecay.ProposeTrackStatus( fStopAndKill ) ;
  fParticleChangeForDecay.ProposeLocalEnergyDeposit(energyDeposit); 
  fParticleChangeForDecay.ProposeGlobalTime( finalGlobalTime );
  // reset NumberOfInteractionLengthLeft
  ClearNumberOfInteractionLengthLeft();

  return &fParticleChangeForDecay ;
} 






