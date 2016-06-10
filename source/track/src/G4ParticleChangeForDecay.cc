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
// $Id: G4ParticleChangeForDecay.cc 93028 2015-09-30 16:09:00Z gcosmo $
//
// 
// --------------------------------------------------------------
//	GEANT 4 class implementation file 
//
//	
//	
// ------------------------------------------------------------
//   Implemented for the new scheme                 23 Mar. 1998  H.Kurahige
//   Remove modification of energy/momentum         20 Jul, 1998  H.Kurashige
// --------------------------------------------------------------

#include "G4ParticleChangeForDecay.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4TrackFastVector.hh"
#include "G4DynamicParticle.hh"
#include "G4ExceptionSeverity.hh"

G4ParticleChangeForDecay::G4ParticleChangeForDecay()
  : G4VParticleChange(), 
    theGlobalTime0(0.),theLocalTime0(0.),theTimeChange(0.)
{
#ifdef G4VERBOSE
  if (verboseLevel>2) {
    G4cout << "G4ParticleChangeForDecay::G4ParticleChangeForDecay() " << G4endl;
  }
#endif
}

G4ParticleChangeForDecay::~G4ParticleChangeForDecay() 
{
#ifdef G4VERBOSE
  if (verboseLevel>2) {
    G4cout << "G4ParticleChangeForDecay::~G4ParticleChangeForDecay() " << G4endl;
  }
#endif
}

G4ParticleChangeForDecay::G4ParticleChangeForDecay(const G4ParticleChangeForDecay &right): G4VParticleChange(right)
{
  theGlobalTime0 = right.theGlobalTime0;
  theLocalTime0 = right.theLocalTime0;
  theTimeChange = right.theTimeChange;
  thePolarizationChange = right.thePolarizationChange;
}


G4ParticleChangeForDecay & G4ParticleChangeForDecay::operator=(const G4ParticleChangeForDecay &right)
{
  if (this != &right){
    if (theNumberOfSecondaries>0) {
#ifdef G4VERBOSE
      if (verboseLevel>0) {
	G4cout << "G4ParticleChangeForDecay: assignment operator Warning  ";
	G4cout << "theListOfSecondaries is not empty ";
       }
#endif
      for (G4int index= 0; index<theNumberOfSecondaries; index++){
	if ( (*theListOfSecondaries)[index] ) delete (*theListOfSecondaries)[index] ;
      }
    }
    delete theListOfSecondaries; 
    theListOfSecondaries =  new G4TrackFastVector();
    theNumberOfSecondaries = right.theNumberOfSecondaries;
    for (G4int index = 0; index<theNumberOfSecondaries; index++){
      G4Track* newTrack =  new G4Track(*((*right.theListOfSecondaries)[index] ));
      theListOfSecondaries->SetElement(index, newTrack);			    }

    theStatusChange = right.theStatusChange;
    theTrueStepLength = right.theTrueStepLength;
    theLocalEnergyDeposit = right.theLocalEnergyDeposit;
    theSteppingControlFlag = right.theSteppingControlFlag;
    
    theGlobalTime0 = right.theGlobalTime0;
    theLocalTime0 = right.theLocalTime0;
    theTimeChange = right.theTimeChange;
    thePolarizationChange = right.thePolarizationChange;
  }
  return *this;
}

G4bool G4ParticleChangeForDecay::operator==(const G4ParticleChangeForDecay &right) const
{
   return ((G4VParticleChange *)this == (G4VParticleChange *) &right);
}

G4bool G4ParticleChangeForDecay::operator!=(const G4ParticleChangeForDecay &right) const
{
   return ((G4VParticleChange *)this != (G4VParticleChange *) &right);
}

//----------------------------------------------------------------
// methods for Initialization
//
void G4ParticleChangeForDecay::Initialize(const G4Track& track)
{
  // use base class's method at first
  G4VParticleChange::Initialize(track);

  const G4DynamicParticle*  pParticle = track.GetDynamicParticle();

  // set TimeChange equal to local time of the parent track
  theTimeChange                = track.GetLocalTime();

  // set initial Local/Global time of the parent track
  theLocalTime0           = track.GetLocalTime();
  theGlobalTime0          = track.GetGlobalTime();

  // set the Polarization equal to those of the parent track
  thePolarizationChange  = pParticle->GetPolarization();
}

//----------------------------------------------------------------
// methods for updating G4Step 
//

G4Step* G4ParticleChangeForDecay::UpdateStepForPostStep(G4Step* pStep)
{ 
  G4StepPoint* pPostStepPoint = pStep->GetPostStepPoint(); 

  if (isParentWeightProposed ){
    pPostStepPoint->SetWeight( theParentWeight );
  }
  // update polarization
  pPostStepPoint->SetPolarization( thePolarizationChange );

  //  Update the G4Step specific attributes 
  return UpdateStepInfo(pStep);
}


G4Step* G4ParticleChangeForDecay::UpdateStepForAtRest(G4Step* pStep)
{ 
  // A physics process always calculates the final state of the particle

  G4StepPoint* pPostStepPoint = pStep->GetPostStepPoint(); 

  // update polarization
  pPostStepPoint->SetPolarization( thePolarizationChange );
 
  // update time
  pPostStepPoint->SetGlobalTime( GetGlobalTime() );
  pPostStepPoint->SetLocalTime(  theTimeChange ); 
  pPostStepPoint->AddProperTime (theTimeChange-theLocalTime0);

#ifdef G4VERBOSE
  G4Track*     aTrack  = pStep->GetTrack();
  if (debugFlag) CheckIt(*aTrack);
#endif

  if (isParentWeightProposed )pPostStepPoint->SetWeight( theParentWeight );

  //  Update the G4Step specific attributes 
  return UpdateStepInfo(pStep);
}

void G4ParticleChangeForDecay::DumpInfo() const
{
// Show header
  G4VParticleChange::DumpInfo();

  G4int oldprc = G4cout.precision(3);
  G4cout << " proposed local Time (ns)     : " 
         << std::setw(20) << theTimeChange/ns << G4endl;
  G4cout << " initial local Time (ns)      : "
	 << std::setw(20) << theLocalTime0/ns << G4endl;
  G4cout << " initial global Time (ns)      : "
	 << std::setw(20) << theGlobalTime0/ns << G4endl;
  G4cout.precision(oldprc);
}

G4bool G4ParticleChangeForDecay::CheckIt(const G4Track& aTrack)
{
  G4bool    exitWithError = false;

  G4double  accuracy;

  // local time should not go back
  G4bool itsOK =true;
  accuracy = -1.0*(theTimeChange - theLocalTime0)/ns;
  if (accuracy > accuracyForWarning) {
    itsOK = false;
    exitWithError = (accuracy > accuracyForException);
#ifdef G4VERBOSE
    G4cout << "  G4ParticleChangeForDecay::CheckIt    : ";
    G4cout << "the local time goes back  !!" 
	   << "  Difference:  " << accuracy  << "[ns] " <<G4endl;
    G4cout << "initial local time "<< theLocalTime0/ns <<   "[ns] " 
	   << "initial global time "<< theGlobalTime0/ns  << "[ns] " <<G4endl;
    G4cout << aTrack.GetDefinition()->GetParticleName()
	   << " E=" << aTrack.GetKineticEnergy()/MeV
	   << " pos=" << aTrack.GetPosition().x()/m
	   << ", " << aTrack.GetPosition().y()/m
	   << ", " << aTrack.GetPosition().z()/m
	   <<G4endl;
#endif
  }

  // dump out information of this particle change
#ifdef G4VERBOSE
  if (!itsOK) DumpInfo();
#endif

  // Exit with error
  if (exitWithError) {
    G4Exception("G4ParticleChangeForDecay::CheckIt",
		"TRACK005",EventMustBeAborted,
		"time was  illegal");
  } 

  // correction
  if (!itsOK) {
    theTimeChange = aTrack.GetLocalTime();
  }

  itsOK = (itsOK) && G4VParticleChange::CheckIt(aTrack);
  return itsOK;
}





