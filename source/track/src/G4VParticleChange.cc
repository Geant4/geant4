// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VParticleChange.cc,v 1.1 1999-01-07 16:14:27 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//	GEANT 4 class implementation file 
//
//	For information related to this code contact:
//	CERN, CN Division, ASD Group
//
// ------------------------------------------------------------
//   Implemented for the new scheme                 23 Mar. 1998  H.Kurahige
// --------------------------------------------------------------

#include "G4VParticleChange.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4TrackFastVector.hh"
#include "G4Mars5GeVMechanism.hh"

G4VParticleChange::G4VParticleChange():
   theNumberOfSecondaries(0),
   theSizeOftheListOfSecondaries(G4TrackFastVectorSize),
   theStatusChange(fAlive),
   theSteppingControlFlag(NormalCondition),     
   theLocalEnergyDeposit(0.0),
   theParentWeight(1.0),
   theEBMechanism(NULL),
   fUseEB(false),
   verboseLevel(1)
{
   theListOfSecondaries = new G4TrackFastVector();
}

G4VParticleChange::G4VParticleChange(G4bool useEB):
   theNumberOfSecondaries(0),
   theSizeOftheListOfSecondaries(G4TrackFastVectorSize),
   theStatusChange(fAlive),
   theSteppingControlFlag(NormalCondition),     
   theLocalEnergyDeposit(0.0),
   theParentWeight(1.0),
   fUseEB(useEB),
   verboseLevel(1)
{
   theListOfSecondaries = new G4TrackFastVector();
   // register  G4EvtBiasMechanism as a default
   theEBMechanism = new G4Mars5GeVMechanism();
}

G4VParticleChange::~G4VParticleChange() {
  // check if tracks still exist in theListOfSecondaries
  if (theNumberOfSecondaries>0) {
#ifdef G4VERBOSE
    if (verboseLevel>0) {
      G4cerr << "G4VParticleChange::~G4VParticleChange() Warning  ";
      G4cerr << "theListOfSecondaries is not empty ";
    }
#endif
    for (G4int index= 0; index<theNumberOfSecondaries; index++){
      if ( (*theListOfSecondaries)[index] ) delete (*theListOfSecondaries)[index] ;
    }
  }
  if (theEBMechanism !=NULL) delete theEBMechanism;
  delete theListOfSecondaries; 
}

// copy and assignment operators are implemented as "shallow copy"
G4VParticleChange::G4VParticleChange(const G4VParticleChange &right)
{
   *this = right;
}


G4VParticleChange & G4VParticleChange::operator=(const G4VParticleChange &right)
{
   if (this != &right)
   {
      theListOfSecondaries = right.theListOfSecondaries;
      theSizeOftheListOfSecondaries = right.theSizeOftheListOfSecondaries;
      theNumberOfSecondaries = right.theNumberOfSecondaries;
      theStatusChange = right.theStatusChange;
      theTrueStepLength = right.theTrueStepLength;
      theLocalEnergyDeposit = right.theLocalEnergyDeposit;
      theSteppingControlFlag = right.theSteppingControlFlag;
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

//----------------------------------------------------------------
// methods for printing messages  
//
 
void G4VParticleChange::DumpInfo() const
{

// Show header
  G4cout.precision(3);
  G4cout << "      -----------------------------------------------" 
       << endl;
  G4cout << "        G4ParticleChange Information  " << setw(20) << endl;
  G4cout << "      -----------------------------------------------" 
       << endl;

  G4cout << "        # of 2ndaries       : " 
       << setw(20) << theNumberOfSecondaries
       << endl;

  if (theNumberOfSecondaries >0) {
    G4cout << "        Pointer to 2ndaries : " 
         << setw(20) << GetSecondary(0)
         << endl;
    G4cout << "        (Showed only 1st one)"
         << endl;
  }
  G4cout << "      -----------------------------------------------" 
       << endl;

  G4cout << "        Energy Deposit (MeV): " 
       << setw(20) << theLocalEnergyDeposit/MeV
       << endl;
  G4cout << "        Track Status        : " 
       << setw(20);
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
       G4cout << endl;
  G4cout << "        True Path Length (mm) : " 
       << setw(20) << theTrueStepLength/mm
       << endl;
  G4cout << "        Stepping Control     : " 
       << setw(20) << theSteppingControlFlag
       << endl;   
  G4cout << "        Event Biasing        : ";
  if (fUseEB) {
    G4cout << setw(20) << theEBMechanism->GetName();
  } else {
    G4cout << " not used ";
  }
  G4cout << endl;      
}







