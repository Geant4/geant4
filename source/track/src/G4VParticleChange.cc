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
// $Id: G4VParticleChange.cc 92776 2015-09-16 06:57:55Z gcosmo $
//
// 
// --------------------------------------------------------------
//	GEANT 4 class implementation file 
//
//
// ------------------------------------------------------------
//   Implemented for the new scheme                 23 Mar. 1998  H.Kurahige
// --------------------------------------------------------------

#include "G4VParticleChange.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4TrackFastVector.hh"
#include "G4ExceptionSeverity.hh"

const G4double G4VParticleChange::accuracyForWarning = 1.0e-9;
const G4double G4VParticleChange::accuracyForException = 0.001;

G4VParticleChange::G4VParticleChange()
  :theListOfSecondaries(0),
   theNumberOfSecondaries(0),
   theSizeOftheListOfSecondaries(G4TrackFastVectorSize),
   theStatusChange(fAlive),
   theSteppingControlFlag(NormalCondition),     
   theLocalEnergyDeposit(0.0),
   theNonIonizingEnergyDeposit(0.0),
   theTrueStepLength(0.0),
   theFirstStepInVolume(false),
   theLastStepInVolume(false),
   theParentWeight(1.0),
   isParentWeightProposed(false),
   fSetSecondaryWeightByProcess(false),
   theParentGlobalTime(0.0),
   verboseLevel(1), 
   debugFlag(false)
{
#ifdef G4VERBOSE
  // activate CHeckIt if in VERBOSE mode
  debugFlag = true;
#endif
  theListOfSecondaries = new G4TrackFastVector();
}

G4VParticleChange::~G4VParticleChange() {
  // check if tracks still exist in theListOfSecondaries
  if (theNumberOfSecondaries>0) {
#ifdef G4VERBOSE
    if (verboseLevel>0) {
      G4cout << "G4VParticleChange::~G4VParticleChange() Warning  ";
      G4cout << "theListOfSecondaries is not empty ";
    }
#endif
    for (G4int index= 0; index<theNumberOfSecondaries; index++){
      delete (*theListOfSecondaries)[index] ;
    }
  }
  delete theListOfSecondaries; 
}

G4VParticleChange::G4VParticleChange(const G4VParticleChange &right)
  :theListOfSecondaries(0),
   theNumberOfSecondaries(0),
   theSizeOftheListOfSecondaries(G4TrackFastVectorSize),
   theStatusChange( right.theStatusChange),
   theSteppingControlFlag(right.theSteppingControlFlag),
   theLocalEnergyDeposit(right.theLocalEnergyDeposit),
   theNonIonizingEnergyDeposit(right.theNonIonizingEnergyDeposit),
   theTrueStepLength(right.theTrueStepLength),
   theFirstStepInVolume( right.theFirstStepInVolume),
   theLastStepInVolume(right.theLastStepInVolume),
   theParentWeight(right.theParentWeight),
   isParentWeightProposed(false),
   fSetSecondaryWeightByProcess(right.fSetSecondaryWeightByProcess),
   theParentGlobalTime(0.0),
   verboseLevel(right.verboseLevel),
   debugFlag(right.debugFlag)
{
  theListOfSecondaries =  new G4TrackFastVector();
  theNumberOfSecondaries = right.theNumberOfSecondaries;
  for (G4int index = 0; index<theNumberOfSecondaries; index++){
    G4Track* newTrack =  new G4Track(*((*right.theListOfSecondaries)[index] ));
    theListOfSecondaries->SetElement(index, newTrack);			    
  }
}


G4VParticleChange & G4VParticleChange::operator=(const G4VParticleChange &right)
{
  if (this != &right){
    if (theNumberOfSecondaries>0) {
#ifdef G4VERBOSE
      if (verboseLevel>0) {
	G4cout << "G4VParticleChange: assignment operator Warning  ";
	G4cout << "theListOfSecondaries is not empty ";
      }
#endif
      for (G4int index = 0; index<theNumberOfSecondaries; index++){
	if ( (*theListOfSecondaries)[index] ) delete ((*theListOfSecondaries)[index]) ;
      }
    }
    delete theListOfSecondaries; 
      
    theListOfSecondaries =  new G4TrackFastVector();
    theNumberOfSecondaries = right.theNumberOfSecondaries;
   for (G4int index = 0; index<theNumberOfSecondaries; index++){
    G4Track* newTrack =  new G4Track(*((*right.theListOfSecondaries)[index] ));
    theListOfSecondaries->SetElement(index, newTrack);			    
  }
    theStatusChange = right.theStatusChange;
    theSteppingControlFlag = right.theSteppingControlFlag;
    theLocalEnergyDeposit = right.theLocalEnergyDeposit;
    theNonIonizingEnergyDeposit = right.theNonIonizingEnergyDeposit;
    theTrueStepLength = right.theTrueStepLength;
    
    theFirstStepInVolume = right.theFirstStepInVolume;
    theLastStepInVolume =  right.theLastStepInVolume;

    theParentWeight = right.theParentWeight;
    isParentWeightProposed = right.isParentWeightProposed;
    fSetSecondaryWeightByProcess = right.fSetSecondaryWeightByProcess;

    theParentGlobalTime = right.theParentGlobalTime;

    verboseLevel = right.verboseLevel;
    debugFlag = right.debugFlag;
    
  }
  return *this;
}

G4bool G4VParticleChange::operator==(const G4VParticleChange &right) const
{
   return (this == (G4VParticleChange *) &right);
}


G4bool G4VParticleChange::operator!=(const G4VParticleChange &right) const
{
   return (this != (G4VParticleChange *) &right);
}

void G4VParticleChange::AddSecondary(G4Track *aTrack)
{
  if (debugFlag) CheckSecondary(*aTrack);

  // add a secondary after size check
  if (theSizeOftheListOfSecondaries > theNumberOfSecondaries) {
    // Set weight of secondary tracks
    if (!fSetSecondaryWeightByProcess) aTrack->SetWeight(theParentWeight);
    theListOfSecondaries->SetElement(theNumberOfSecondaries, aTrack);
    theNumberOfSecondaries++;
  } else {
    delete aTrack;

#ifdef G4VERBOSE
    if (verboseLevel>0) {
      G4cout << "G4VParticleChange::AddSecondary() Warning  ";
      G4cout << "theListOfSecondaries is full !! " << G4endl;
      G4cout << " The track is deleted " << G4endl;
    }
#endif
    G4Exception("G4VParticleChange::AddSecondary",
                "TRACK101", JustWarning,
                "Secondary Bug is full. The track is deleted"); 
  }
}


 
// Virtual methods for updating G4Step 
//

G4Step* G4VParticleChange::UpdateStepInfo(G4Step* pStep)
{
  // Update the G4Step specific attributes
  pStep->SetStepLength( theTrueStepLength );
  pStep->AddTotalEnergyDeposit( theLocalEnergyDeposit );
  pStep->AddNonIonizingEnergyDeposit( theNonIonizingEnergyDeposit );
  pStep->SetControlFlag( theSteppingControlFlag );

  if (theFirstStepInVolume) {pStep->SetFirstStepFlag();}
  else                      {pStep->ClearFirstStepFlag();} 
  if (theLastStepInVolume)  {pStep->SetLastStepFlag();}
  else                      {pStep->ClearLastStepFlag();} 

  return pStep;
}


G4Step* G4VParticleChange::UpdateStepForAtRest(G4Step* Step)
{ 
  if (isParentWeightProposed ){
    Step->GetPostStepPoint()->SetWeight( theParentWeight );
  }
  return UpdateStepInfo(Step);
}


G4Step* G4VParticleChange::UpdateStepForAlongStep(G4Step* Step)
{
  if (isParentWeightProposed ){ 
    G4double initialWeight = Step->GetPreStepPoint()->GetWeight();
    G4double currentWeight = Step->GetPostStepPoint()->GetWeight();
    G4double finalWeight   = (theParentWeight/initialWeight)*currentWeight;
    Step->GetPostStepPoint()->SetWeight( finalWeight );
  }   
  return UpdateStepInfo(Step);
}

G4Step* G4VParticleChange::UpdateStepForPostStep(G4Step* Step)
{
  if (isParentWeightProposed ){ 
    Step->GetPostStepPoint()->SetWeight( theParentWeight );
  }
  return UpdateStepInfo(Step);
}


//----------------------------------------------------------------
// methods for printing messages  
//
 
void G4VParticleChange::DumpInfo() const
{

// Show header
  G4int olprc = G4cout.precision(3);
  G4cout << "      -----------------------------------------------" 
       << G4endl;
  G4cout << "        G4ParticleChange Information  " << std::setw(20) << G4endl;
  G4cout << "      -----------------------------------------------" 
       << G4endl;

  G4cout << "        # of 2ndaries       : " 
       << std::setw(20) << theNumberOfSecondaries
       << G4endl;

  if (theNumberOfSecondaries >0) {
    G4cout << "        Pointer to 2ndaries : " 
         << std::setw(20) << GetSecondary(0)
         << G4endl;
    G4cout << "        (Showed only 1st one)"
         << G4endl;
  }
  G4cout << "      -----------------------------------------------" 
       << G4endl;

  G4cout << "        Energy Deposit (MeV): " 
       << std::setw(20) << theLocalEnergyDeposit/MeV
       << G4endl;

  G4cout << "        Non-ionizing Energy Deposit (MeV): " 
       << std::setw(20) << theNonIonizingEnergyDeposit/MeV
       << G4endl;

  G4cout << "        Track Status        : " 
       << std::setw(20);
       if( theStatusChange == fAlive ){
         G4cout << " Alive";
       } else if( theStatusChange == fStopButAlive ){
           G4cout << " StopButAlive";
       } else if( theStatusChange == fStopAndKill ){
           G4cout << " StopAndKill";
       } else if( theStatusChange  == fKillTrackAndSecondaries ){
           G4cout << " KillTrackAndSecondaries";
       } else if( theStatusChange  == fSuspend ){
           G4cout << " Suspend";
       } else if( theStatusChange == fPostponeToNextEvent ){
           G4cout << " PostponeToNextEvent";
       }
       G4cout << G4endl;
  G4cout << "        True Path Length (mm) : " 
       << std::setw(20) << theTrueStepLength/mm
       << G4endl;
  G4cout << "        Stepping Control      : " 
       << std::setw(20) << theSteppingControlFlag
       << G4endl;   
  if (theFirstStepInVolume) {
    G4cout << "    First Step In the voulme  : "  << G4endl;
  }
  if (theLastStepInVolume) {
    G4cout << "    Last Step In the voulme  : "  << G4endl;
  }
  G4cout.precision(olprc);
}

G4bool G4VParticleChange::CheckIt(const G4Track&
#ifdef G4VERBOSE
                                         aTrack
#endif
)
{

  G4bool    exitWithError = false;
  G4double  accuracy;
  static G4ThreadLocal G4int nError = 0;
#ifdef G4VERBOSE
  const  G4int maxError = 30;
#endif

  // Energy deposit should not be negative
  G4bool itsOKforEnergy = true;
  accuracy = -1.0*theLocalEnergyDeposit/MeV;
  if (accuracy > accuracyForWarning) {
    itsOKforEnergy = false;
    nError += 1;
    exitWithError =  (accuracy > accuracyForException);
#ifdef G4VERBOSE
    if (nError < maxError) {
      G4cout << "  G4VParticleChange::CheckIt    : ";
      G4cout << "the energy deposit  is negative  !!" 
	     << "  Difference:  " << accuracy  << "[MeV] " <<G4endl;
      G4cout << aTrack.GetDefinition()->GetParticleName()
	     << " E=" << aTrack.GetKineticEnergy()/MeV
	     << " pos=" << aTrack.GetPosition().x()/m
	     << ", " << aTrack.GetPosition().y()/m
	     << ", " << aTrack.GetPosition().z()/m
	     <<G4endl;
    }
#endif
  }
 
  // true path length should not be negative
  G4bool itsOKforStepLength = true;
  accuracy = -1.0*theTrueStepLength/mm;
  if (accuracy > accuracyForWarning) {
    itsOKforStepLength = false;
    nError += 1;
    exitWithError =  (accuracy > accuracyForException);
#ifdef G4VERBOSE
    if (nError < maxError) {
      G4cout << "  G4VParticleChange::CheckIt    : ";
      G4cout << "the true step length is negative  !!"
	     << "  Difference:  " << accuracy  << "[MeV] " <<G4endl;
      G4cout << aTrack.GetDefinition()->GetParticleName()
	     << " E=" << aTrack.GetKineticEnergy()/MeV
	     << " pos=" << aTrack.GetPosition().x()/m
	     << ", " << aTrack.GetPosition().y()/m
	     << ", " << aTrack.GetPosition().z()/m
	     <<G4endl;
    }
#endif

  }
#ifdef G4VERBOSE
  if (!itsOKforStepLength || !itsOKforEnergy ){
  // dump out information of this particle change
    DumpInfo();
  }
#endif

  // Exit with error
  if (exitWithError) {
    G4Exception("G4VParticleChange::CheckIt",
		"TRACK001", EventMustBeAborted,
		"Step length and/or energy deposit was illegal");
  }

  // correction 
  if ( !itsOKforStepLength ) {
    theTrueStepLength =  (1.e-12)*mm;
  } 
  if ( !itsOKforEnergy ) {
    theLocalEnergyDeposit = 0.0;
  }
  return (itsOKforStepLength && itsOKforEnergy );
}

G4bool G4VParticleChange::CheckSecondary(G4Track& aTrack)
{
  G4bool    exitWithError = false;
  G4double  accuracy;
  static G4ThreadLocal G4int nError = 0;
#ifdef G4VERBOSE
  const  G4int maxError = 30;
#endif

  // MomentumDirection should be unit vector
  G4bool itsOKforMomentum = true;  
  if (aTrack.GetKineticEnergy()>0.){
    accuracy = std::fabs((aTrack.GetMomentumDirection()).mag2()-1.0);
    if (accuracy > accuracyForWarning) {
      itsOKforMomentum = false;
      nError += 1;
      exitWithError = exitWithError || (accuracy > accuracyForException);
#ifdef G4VERBOSE
      if (nError < maxError) {
	G4cout << " G4VParticleChange::CheckSecondary  :   ";
	G4cout << "the Momentum direction is not unit vector !! " 
	       << "  Difference:  " << accuracy << G4endl;
	G4cout << aTrack.GetDefinition()->GetParticleName()
	       << " E=" << aTrack.GetKineticEnergy()/MeV
	       << " pos=" << aTrack.GetPosition().x()/m
	       << ", " << aTrack.GetPosition().y()/m
	       << ", " << aTrack.GetPosition().z()/m
	       <<G4endl;
      }
#endif
    }
  }
  
  // Kinetic Energy should not be negative
  G4bool itsOKforEnergy = true;
  accuracy = -1.0*(aTrack.GetKineticEnergy())/MeV;
  if (accuracy > accuracyForWarning) {
    itsOKforEnergy = false;
    nError += 1;
    exitWithError = exitWithError ||  (accuracy > accuracyForException);
#ifdef G4VERBOSE
    if (nError < maxError) {
      G4cout << " G4VParticleChange::CheckSecondary  :   ";
      G4cout << "the kinetic energy is negative  !!" 
	     << "  Difference:  " << accuracy  << "[MeV] " <<G4endl;
      G4cout << " G4VParticleChange::CheckSecondary  :   ";
      G4cout << "the global time of secondary is earlier than the parent  !!" 
	     << "  Difference:  " << accuracy  << "[ns] " <<G4endl;
      G4cout << aTrack.GetDefinition()->GetParticleName()
	     << " E=" << aTrack.GetKineticEnergy()/MeV
	     << " pos=" << aTrack.GetPosition().x()/m
	     << ", " << aTrack.GetPosition().y()/m
	     << ", " << aTrack.GetPosition().z()/m
	   <<G4endl;
    }
#endif
  }
  // Check timing of secondaries
  G4bool itsOKforTiming = true;

  accuracy = (theParentGlobalTime - aTrack.GetGlobalTime())/ns;
  if (accuracy > accuracyForWarning){
    itsOKforTiming = false;
    nError += 1;
    exitWithError = (accuracy > accuracyForException);
#ifdef G4VERBOSE
    if (nError < maxError) {
      G4cout << " G4VParticleChange::CheckSecondary  :   ";
      G4cout << "the global time of secondary goes back comapared to the parent  !!" 
	     << "  Difference:  " << accuracy  << "[ns] " <<G4endl;
      G4cout << aTrack.GetDefinition()->GetParticleName()
	     << " E=" << aTrack.GetKineticEnergy()/MeV
	     << " pos=" << aTrack.GetPosition().x()/m
	     << ", " << aTrack.GetPosition().y()/m
	       << ", " << aTrack.GetPosition().z()/m
	     << " time=" << aTrack.GetGlobalTime()/ns
	     << " parent time=" << theParentGlobalTime/ns
	     <<G4endl;
    }
#endif
  }

  // Exit with error
  if (exitWithError) {
    G4Exception("G4VParticleChange::CheckSecondary",
		"TRACK001", EventMustBeAborted,
		"Secondary with illegal energy/momentum ");
  }

  G4bool itsOK = itsOKforMomentum && itsOKforEnergy && itsOKforTiming;

  //correction
  if (!itsOKforMomentum) {
    G4double vmag = (aTrack.GetMomentumDirection()).mag();
    aTrack.SetMomentumDirection((1./vmag)*aTrack.GetMomentumDirection());
  }
  if (!itsOKforEnergy) {
    aTrack.SetKineticEnergy(0.0);
  }
 
  if (!itsOK) {
    this->DumpInfo();
    
  }


  return itsOK;
}


G4double G4VParticleChange::GetAccuracyForWarning() const
{
  return accuracyForWarning;
}

G4double G4VParticleChange::GetAccuracyForException() const
{
  return accuracyForException;
}


///////////////////////////////////////////////////////////
//Obsolete methods for parent weight
/////////////////////
void  G4VParticleChange::SetParentWeightByProcess(G4bool )
{
}


G4bool   G4VParticleChange::IsParentWeightSetByProcess() const
{
  return  true;
}
