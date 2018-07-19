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
// $Id: G4UnknownDecay.cc 105727 2017-08-16 12:47:05Z gcosmo $
//
// 
// --------------------------------------------------------------
//      GEANT 4 class implementation file
//
// ------------------------------------------------------------
//

#include "G4UnknownDecay.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4DynamicParticle.hh"
#include "G4DecayProducts.hh"
#include "G4PhysicsLogVector.hh"
#include "G4ParticleChangeForDecay.hh"
#include "G4DecayProcessType.hh"

// constructor
G4UnknownDecay::G4UnknownDecay(const G4String& processName)
                               :G4VDiscreteProcess(processName, fDecay),
				verboseLevel(1),
                                HighestValue(20.0)
{
  // set Process Sub Type
  SetProcessSubType(static_cast<int>(DECAY_Unknown));

#ifdef G4VERBOSE
  if (GetVerboseLevel()>1) {
    G4cout << "G4UnknownDecay  constructor " << "  Name:" << processName << G4endl;
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
  G4bool isPreAssigned = (o_products != nullptr);   
  G4DecayProducts* products = nullptr;

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
  G4double   ParentMass    = aParticle->GetMass();
  if (ParentEnergy < ParentMass) {
    ParentEnergy = ParentMass;
#ifdef G4VERBOSE
    if (GetVerboseLevel()>1) {
      G4cout << "G4UnknownDecay::DoIt  : Total Energy is less than its mass" << G4endl;
      G4cout << " Particle: " << aParticle->GetDefinition()->GetParticleName();
      G4cout << " Energy:"    << ParentEnergy/MeV << "[MeV]";
      G4cout << " Mass:"    << ParentMass/MeV << "[MeV]";
      G4cout << G4endl;
    }
#endif
  }
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
    G4cout << "G4UnknownDecay::DoIt  : Decay vertex :";
    G4cout << " Time: " << finalGlobalTime/ns << "[ns]";
    G4cout << " X:" << (aTrack.GetPosition()).x() /cm << "[cm]";
    G4cout << " Y:" << (aTrack.GetPosition()).y() /cm << "[cm]";
    G4cout << " Z:" << (aTrack.GetPosition()).z() /cm << "[cm]";
    G4cout << G4endl;
    G4cout << "G4UnknownDecay::DoIt  : decay products in Lab. Frame" << G4endl;
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
  fParticleChangeForDecay.ProposeGlobalTime( finalGlobalTime );
  // reset NumberOfInteractionLengthLeft
  ClearNumberOfInteractionLengthLeft();

  return &fParticleChangeForDecay ;
} 

void G4UnknownDecay::ProcessDescription(std::ostream& outFile) const
{
  outFile << GetProcessName()
	  << ": Decay of 'unknown' particles. \n"
	  << "kinematics of daughters are dertermined "
	  << "by PreAssignedDecayProducts. \n";
}





