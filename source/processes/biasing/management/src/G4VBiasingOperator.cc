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
#include "G4VBiasingOperator.hh"
#include "G4VBiasingOperation.hh"
#include "G4VParticleChange.hh"


G4MapCache< const G4LogicalVolume*, G4VBiasingOperator* > G4VBiasingOperator::fLogicalToSetupMap;
G4VectorCache< G4VBiasingOperator* >                      G4VBiasingOperator::fOperators;
G4Cache< G4BiasingOperatorStateNotifier* >                G4VBiasingOperator::fStateNotifier(0);


G4VBiasingOperator::G4VBiasingOperator(G4String name)
  : fName                                      ( name     ),
    fOccurenceBiasingOperation                 ( nullptr  ),
    fFinalStateBiasingOperation                ( nullptr  ),
    fNonPhysicsBiasingOperation                ( nullptr  ),
    fPreviousProposedOccurenceBiasingOperation ( nullptr  ),
    fPreviousProposedFinalStateBiasingOperation( nullptr  ),
    fPreviousProposedNonPhysicsBiasingOperation( nullptr  ),
    fPreviousAppliedOccurenceBiasingOperation  ( nullptr  ),
    fPreviousAppliedFinalStateBiasingOperation ( nullptr  ),
    fPreviousAppliedNonPhysicsBiasingOperation ( nullptr  ),
    fPreviousBiasingAppliedCase                ( BAC_None )
{
  fOperators.Push_back(this);

  if ( fStateNotifier.Get() == 0 ) fStateNotifier.Put( new G4BiasingOperatorStateNotifier() );

}

G4VBiasingOperator::~G4VBiasingOperator()
{
}

void G4VBiasingOperator::AttachTo(const G4LogicalVolume* logical)
{
  G4MapCache< const G4LogicalVolume*, G4VBiasingOperator* >::iterator it;
  it = fLogicalToSetupMap.Find(logical);
  if ( it == fLogicalToSetupMap.End() ) fLogicalToSetupMap[logical] = this;
  else if ( (*it).second != this )
    {
      G4ExceptionDescription ed;
      ed << "Biasing operator `" << GetName() 
	 << "' can not be attached to Logical volume `"
	 << logical->GetName() << "' which is already used by another operator !" << G4endl;
      G4Exception("G4VBiasingOperator::AttachTo(...)",
		  "BIAS.MNG.01",
		  JustWarning,
		  ed);
    }
}


G4VBiasingOperator* G4VBiasingOperator::GetBiasingOperator(const G4LogicalVolume* logical)
{
  G4MapCache< const G4LogicalVolume*, G4VBiasingOperator* >::iterator it;
  it = fLogicalToSetupMap.Find(logical);
  if ( it == fLogicalToSetupMap.End() ) return nullptr;
  else return (*it).second;
}

G4VBiasingOperation* G4VBiasingOperator::GetProposedOccurenceBiasingOperation(const G4Track* track, const G4BiasingProcessInterface* callingProcess)
{
  fOccurenceBiasingOperation = ProposeOccurenceBiasingOperation(track, callingProcess);
  return fOccurenceBiasingOperation;
}

G4VBiasingOperation* G4VBiasingOperator::GetProposedFinalStateBiasingOperation(const G4Track* track, const G4BiasingProcessInterface* callingProcess)
{
  fFinalStateBiasingOperation = ProposeFinalStateBiasingOperation(track, callingProcess); 
  return fFinalStateBiasingOperation;
}

G4VBiasingOperation* G4VBiasingOperator::GetProposedNonPhysicsBiasingOperation(const G4Track* track, const G4BiasingProcessInterface* callingProcess)
{
  fNonPhysicsBiasingOperation = ProposeNonPhysicsBiasingOperation(track, callingProcess);
  return fNonPhysicsBiasingOperation;
}

void G4VBiasingOperator::ReportOperationApplied( const G4BiasingProcessInterface*  callingProcess,
						 G4BiasingAppliedCase                 biasingCase,
						 G4VBiasingOperation*             operationApplied,
						 const G4VParticleChange*   particleChangeProduced )
{
  fPreviousBiasingAppliedCase = biasingCase;
  fPreviousAppliedOccurenceBiasingOperation  = nullptr;
  fPreviousAppliedFinalStateBiasingOperation = nullptr;
  fPreviousAppliedNonPhysicsBiasingOperation = nullptr;
  switch ( biasingCase )
    {
    case BAC_None:
      break;
    case BAC_NonPhysics:
      fPreviousAppliedNonPhysicsBiasingOperation = operationApplied ;
      break;
    case BAC_FinalState:
      fPreviousAppliedFinalStateBiasingOperation = operationApplied;
      break;
    case BAC_Occurence:
      G4Exception("G4VBiasingOperator::ReportOperationApplied(...)",
      		  "BIAS.MNG.02",
      		  JustWarning,
		  "Internal logic error, please report !");
      break;
    default:
      G4Exception("G4VBiasingOperator::ReportOperationApplied(...)",
		  "BIAS.MNG.03",
		  JustWarning,
		  "Internal logic error, please report !");
    }
  OperationApplied( callingProcess, biasingCase, operationApplied, particleChangeProduced );
}

void G4VBiasingOperator::ReportOperationApplied( const G4BiasingProcessInterface*               callingProcess,
						 G4BiasingAppliedCase                              biasingCase,
						 G4VBiasingOperation*                occurenceOperationApplied,
						 G4double                        weightForOccurenceInteraction,
						 G4VBiasingOperation*               finalStateOperationApplied,
						 const G4VParticleChange*               particleChangeProduced )
{
  fPreviousBiasingAppliedCase = biasingCase;
  fPreviousAppliedOccurenceBiasingOperation  =  occurenceOperationApplied;
  fPreviousAppliedFinalStateBiasingOperation = finalStateOperationApplied;
  OperationApplied( callingProcess, biasingCase, occurenceOperationApplied, weightForOccurenceInteraction, finalStateOperationApplied, particleChangeProduced );
}


void G4VBiasingOperator::ExitingBiasing( const G4Track* track, const G4BiasingProcessInterface* callingProcess )
{
  ExitBiasing( track, callingProcess );
  
  // -- reset all data members:
  fOccurenceBiasingOperation                  = nullptr ;
  fFinalStateBiasingOperation                 = nullptr ;
  fNonPhysicsBiasingOperation                 = nullptr ;
  fPreviousProposedOccurenceBiasingOperation  = nullptr ;
  fPreviousProposedFinalStateBiasingOperation = nullptr ;
  fPreviousProposedNonPhysicsBiasingOperation = nullptr ;
  fPreviousAppliedOccurenceBiasingOperation   = nullptr ;
  fPreviousAppliedFinalStateBiasingOperation  = nullptr ;
  fPreviousAppliedNonPhysicsBiasingOperation  = nullptr ;
  fPreviousBiasingAppliedCase                 = BAC_None ;
}


// -- dummy empty implementations to allow letting arguments visible in the .hh
// -- but avoiding annoying warning messages about unused variables
// -- methods to inform operator that its biasing control is over:
void G4VBiasingOperator::ExitBiasing( const G4Track*, const G4BiasingProcessInterface*)
{}
void G4VBiasingOperator::OperationApplied( const G4BiasingProcessInterface*, G4BiasingAppliedCase,
					   G4VBiasingOperation*, const G4VParticleChange* )
{
}
void G4VBiasingOperator::OperationApplied( const G4BiasingProcessInterface*, G4BiasingAppliedCase,
					   G4VBiasingOperation*,  G4double,
					   G4VBiasingOperation*, const G4VParticleChange* )
{
}


// ----------------------------------------------------------------------------
// -- state machine to get biasing operators messaged at the beginning of runs:
// ----------------------------------------------------------------------------

G4BiasingOperatorStateNotifier::G4BiasingOperatorStateNotifier()
: G4VStateDependent()
{
  fPreviousState =  G4State_PreInit;
}

G4BiasingOperatorStateNotifier::~G4BiasingOperatorStateNotifier()
{}

G4bool G4BiasingOperatorStateNotifier::Notify( G4ApplicationState requestedState )
{
  if ( ( fPreviousState == G4State_Idle ) && ( requestedState == G4State_GeomClosed ) )
    {
      for ( G4int i = 0; i < (G4int)G4VBiasingOperator::fOperators.Size(); ++i )
        G4VBiasingOperator::fOperators[i]->StartRun();
    }

  fPreviousState = requestedState;
  
  return true;
}
