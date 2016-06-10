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
// $Id: G4DNASecondOrderReaction.cc 94218 2015-11-09 08:24:48Z gcosmo $
//
#include "G4DNASecondOrderReaction.hh"

#include <G4VScheduler.hh>
#include "G4SystemOfUnits.hh"
#include "G4Molecule.hh"
#include "G4DNAMolecularMaterial.hh"
#include "G4MolecularConfiguration.hh"
#include "G4DNADamages.hh"
#include "G4UnitsTable.hh"
#include "G4TrackingInformation.hh"

#ifndef State
#define State(theXInfo) (GetState<SecondOrderReactionState>()->theXInfo)
#endif

void G4DNASecondOrderReaction::Create()
{
  pParticleChange = &fParticleChange;
  enableAtRestDoIt    = false;
  enableAlongStepDoIt = false;
  enablePostStepDoIt  = true;

  SetProcessSubType(60);

  G4VITProcess::SetInstantiateProcessState(false);
  // meaning G4DNASecondOrderReaction contains a class inheriting from G4ProcessState

  fIsInitialized = false;
  fpMolecularConfiguration = 0;
  fpMaterial = 0;
  fReactionRate = -1.;
  fConcentration = -1.;
  fMolarMassOfMaterial = -1.;
  fProposesTimeStep = true;
  fReturnedValue = DBL_MAX;
  fpMoleculeDensity = 0;

  verboseLevel = 0;
}

G4DNASecondOrderReaction::G4DNASecondOrderReaction(const G4String &aName, G4ProcessType type) :
        G4VITProcess(aName,type)
{
  Create();
}

G4DNASecondOrderReaction::G4DNASecondOrderReaction(const G4DNASecondOrderReaction& rhs):
        G4VITProcess(rhs)
{
  Create();
}

G4DNASecondOrderReaction::~G4DNASecondOrderReaction()
{
  ;
}
G4DNASecondOrderReaction& G4DNASecondOrderReaction::operator=(const G4DNASecondOrderReaction& rhs)
{
  if (this == &rhs) return *this; // handle self assignment

  //assignment operator
  return *this;
}

G4DNASecondOrderReaction::SecondOrderReactionState::SecondOrderReactionState() : G4ProcessState()
{
  fPreviousTimeAtPreStepPoint = -1;
  fIsInGoodMaterial = false;
}

void G4DNASecondOrderReaction::BuildPhysicsTable(const G4ParticleDefinition&)
{
  fpMoleculeDensity = G4DNAMolecularMaterial::Instance()->GetNumMolPerVolTableFor(fpMaterial);
  fMolarMassOfMaterial = fpMaterial->GetMassOfMolecule()*CLHEP::Avogadro*1e3;
  fIsInitialized = true;
}

void
G4DNASecondOrderReaction::StartTracking(G4Track* track)
{
  G4VProcess::StartTracking(track);
  G4VITProcess::fpState.reset(new SecondOrderReactionState());
  G4VITProcess::StartTracking(track);
}

void
G4DNASecondOrderReaction::SetReaction(const G4MolecularConfiguration* molConf,
                                      const G4Material* mat, double reactionRate)
{
  if(fIsInitialized)
  {
    G4ExceptionDescription exceptionDescription ;
    exceptionDescription << "G4DNASecondOrderReaction was already initialised. ";
    exceptionDescription << "You cannot set a reaction after initialisation.";
    G4Exception("G4DNASecondOrderReaction::SetReaction","G4DNASecondOrderReaction001",
                FatalErrorInArgument,exceptionDescription);
  }
  fpMolecularConfiguration = molConf;
  fpMaterial = mat;
  fReactionRate = reactionRate;
}

G4double G4DNASecondOrderReaction::PostStepGetPhysicalInteractionLength(const G4Track& track,
                                                                        G4double   /*previousStepSize*/,
                                                                        G4ForceCondition* pForceCond)
{
  //    G4cout << "G4DNASecondOrderReaction::PostStepGetPhysicalInteractionLength" << G4endl;
  //    G4cout << "For reaction : " << fpMaterial->GetName() << " + " << fpMolecularConfiguration->GetName() << G4endl;

  //_______________________________________________________________________
  // Check whether the track is in the good material (maybe composite material)
  const G4Material* material = track.GetMaterial();

  G4Molecule* mol = GetMolecule(track);
  if(!mol) return DBL_MAX;
  if(mol->GetMolecularConfiguration() != fpMolecularConfiguration)
  {
    //        G4cout <<"mol->GetMolecularConfiguration() != fpMolecularConfiguration" << G4endl;
    return DBL_MAX;
  }

  G4double molDensity = (*fpMoleculeDensity)[material->GetIndex()];

  if(molDensity == 0.0) // ie : not found
  {
    if(State(fIsInGoodMaterial))
    {
      ResetNumberOfInteractionLengthLeft();
      //State(fPreviousTimeAtPreStepPoint) = -1;
      State(fIsInGoodMaterial) = false;
    }

    // G4cout << " Material " << fpMaterial->GetName() << " not found "
    //        <<" | name of current material : " << material->GetName()
    //        << G4endl;

    return DBL_MAX; // Becareful return here !!
  }

  //    G4cout << " Va calculer le temps d'interaction " << G4endl;

  State(fIsInGoodMaterial) = true;

  //  fConcentration = molDensity/fMolarMassOfMaterial;
  fConcentration = molDensity/CLHEP::Avogadro;
  // G4cout << "Concentration : " << fConcentration / (g/mole)<< G4endl;

  //_______________________________________________________________________
  // Either initialize the lapse of time left
  // meaning
  // => the track enters for the first time in the material
  // or substract the previous time step to the previously calculated lapse
  // of time left
  // meaning
  // => the track has not left this material since the previous call
  G4double previousTimeStep(-1.);

  if(State(fPreviousTimeAtPreStepPoint) != -1)
  {
    previousTimeStep = track.GetGlobalTime() -
        State(fPreviousTimeAtPreStepPoint) ;
  }

  State(fPreviousTimeAtPreStepPoint) = track.GetGlobalTime();

  // condition is set to "Not Forced"
  *pForceCond = NotForced;

  if (
      (previousTimeStep < 0.0) ||
      (fpState->theNumberOfInteractionLengthLeft<=0.0)) {
    // beggining of tracking (or just after DoIt of this process)
    ResetNumberOfInteractionLengthLeft();
  } else if ( previousTimeStep  > 0.0) {
    // get mean free path
    // subtract NumberOfInteractionLengthLeft
    SubtractNumberOfInteractionLengthLeft( previousTimeStep );
  } else {
    // zero time step
    //  Force trigerring the process
    //*pForceCond = Forced;
  }

  fpState->currentInteractionLength = 1/(fReactionRate*fConcentration);

  //    G4cout << "fpState->currentInteractionLength = "
  //        <<  fpState->currentInteractionLength << G4endl;

  G4double value;
  if (fpState->currentInteractionLength <DBL_MAX) {
    value = fpState->theNumberOfInteractionLengthLeft
        * (fpState->currentInteractionLength);
  } else {
    value = DBL_MAX;
  }
#ifdef G4VERBOSE
  if (verboseLevel>2) {
    G4cout << "G4VITRestDiscreteProcess::PostStepGetPhysicalInteractionLength ";
    G4cout << "[ " << GetProcessName() << "]" <<G4endl;
    track.GetDynamicParticle()->DumpInfo();
    G4cout << " in Material  " <<  track.GetMaterial()->GetName() <<G4endl;
    G4cout << "InteractionLength= " << value/cm <<"[cm] " <<G4endl;
  }
#endif

//  G4cout << "currentInteractionLength : " << fpState->currentInteractionLength << G4endl;
//  G4cout << "Material : " << fpMaterial->GetName()
//      << "ID: " << track.GetTrackID()
//      << " Returned time : " << G4BestUnit(value,"Time") << G4endl;

  if(value < fReturnedValue)
    fReturnedValue  = value;

  return value*-1;
  // multiple by -1 to indicate to the tracking system that we are returning a time
}

G4VParticleChange* G4DNASecondOrderReaction::PostStepDoIt(const G4Track& track,const G4Step& /*step*/)
{
  G4Molecule* molecule = GetMolecule(track);
#ifdef G4VERBOSE
  if(verboseLevel > 1)
  {
    G4cout << "___________" << G4endl;
    G4cout << ">>> Beginning of G4DNASecondOrderReaction verbose" << G4endl;
    G4cout << ">>> Returned value : " << G4BestUnit(fReturnedValue,"Time") << G4endl;
    G4cout << ">>> Time Step : " << G4BestUnit(G4VScheduler::Instance()->GetTimeStep(),"Time") << G4endl;
    G4cout << ">>> Reaction : " << molecule->GetName() << " + " << fpMaterial->GetName() << G4endl;
    G4cout << ">>> End of G4DNASecondOrderReaction verbose <<<" << G4endl;
  }
#endif
  fReturnedValue  = DBL_MAX;
  fParticleChange.Initialize(track);
  fParticleChange.ProposeTrackStatus(fStopAndKill);
  G4DNADamages::Instance()->AddIndirectDamage(fpMaterial->GetName(),molecule,track.GetPosition(),track.GetGlobalTime());
  State(fPreviousTimeAtPreStepPoint) = -1;
  return &fParticleChange;
}

