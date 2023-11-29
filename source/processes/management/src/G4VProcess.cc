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
// G4VProcess class implementation
//
// Authors:
// - 2 December 1995, G.Cosmo - First implementation, based on object model
// - 18 December 1996, H.Kurashige - New Physics scheme
// --------------------------------------------------------------------

#include "G4VProcess.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4ProcessTable.hh"
#include "G4PhysicsTable.hh"
#include "G4MaterialTable.hh"
#include "G4ElementTable.hh"
#include "G4ElementVector.hh"
#include "G4Log.hh"

// --------------------------------------------------------------------
G4VProcess::G4VProcess(const G4String& aName, G4ProcessType aType )
  : theProcessName(aName), theProcessType(aType)
{
  pParticleChange = &aParticleChange;
  fProcessTable = G4ProcessTable::GetProcessTable();
  fProcessTable->RegisterProcess(this);
}

// --------------------------------------------------------------------
G4VProcess::G4VProcess()
{
}

// --------------------------------------------------------------------
G4VProcess::~G4VProcess()
{
  fProcessTable->DeRegisterProcess(this);
}

// --------------------------------------------------------------------
G4VProcess::G4VProcess(const G4VProcess& right)
  : theProcessName(right.theProcessName),
    theProcessType(right.theProcessType),
    theProcessSubType(right.theProcessSubType),
    verboseLevel(right.verboseLevel),
    enableAtRestDoIt(right.enableAtRestDoIt),
    enableAlongStepDoIt(right.enableAlongStepDoIt),
    enablePostStepDoIt(right.enablePostStepDoIt),
    masterProcessShadow(right.masterProcessShadow),
    fProcessTable(right.fProcessTable)
{
}

// --------------------------------------------------------------------
void G4VProcess::ResetNumberOfInteractionLengthLeft()
{
  theNumberOfInteractionLengthLeft =  -1.*G4Log( G4UniformRand() );
  theInitialNumberOfInteractionLength = theNumberOfInteractionLengthLeft; 
}

// --------------------------------------------------------------------
void G4VProcess::StartTracking(G4Track*)
{
  currentInteractionLength = -1.0;
  theNumberOfInteractionLengthLeft = -1.0;
  theInitialNumberOfInteractionLength = -1.0;
#ifdef G4VERBOSE
  if (verboseLevel>2)
  {
    G4cout << "G4VProcess::StartTracking() - [" << theProcessName << "]"
           << G4endl;
  }
#endif
}

// --------------------------------------------------------------------
void G4VProcess::EndTracking()
{
#ifdef G4VERBOSE
  if (verboseLevel>2)
  {
    G4cout << "G4VProcess::EndTracking() - [" << theProcessName << "]"
           << G4endl;
  }
#endif
  theNumberOfInteractionLengthLeft = -1.0;
  currentInteractionLength = -1.0;
  theInitialNumberOfInteractionLength=-1.0;
}

// --------------------------------------------------------------------
namespace
{
  static const G4String typeNotDefined = "NotDefined";
  static const G4String typeTransportation = "Transportation";
  static const G4String typeElectromagnetic = "Electromagnetic";
  static const G4String typeOptical = "Optical";
  static const G4String typeHadronic = "Hadronic";
  static const G4String typePhotolepton_hadron = "Photolepton_hadron";
  static const G4String typeDecay = "Decay";
  static const G4String typeGeneral = "General";
  static const G4String typeParameterisation = "Parameterisation";
  static const G4String typeUserDefined = "UserDefined";
  static const G4String typePhonon = "Phonon";
  static const G4String noType = "------";
}

// --------------------------------------------------------------------
const G4String& G4VProcess::GetProcessTypeName(G4ProcessType aType ) 
{
  switch (aType)
  {
    case fNotDefined:         return typeNotDefined; break;
    case fTransportation:     return typeTransportation; break;
    case fElectromagnetic:    return typeElectromagnetic; break;
    case fOptical:            return typeOptical; break;
    case fHadronic:           return typeHadronic; break;
    case fPhotolepton_hadron: return typePhotolepton_hadron; break;
    case fDecay:              return typeDecay; break;
    case fGeneral:            return typeGeneral; break;
    case fParameterisation:   return typeParameterisation; break;
    case fUserDefined:        return typeUserDefined; break;
    case fPhonon:             return typePhonon; break;
    default: ;
  }
  return noType;  
}

// --------------------------------------------------------------------
const G4VProcess* G4VProcess::GetCreatorProcess() const
{
  return this;
}

// --------------------------------------------------------------------
G4bool G4VProcess::operator==(const G4VProcess& right) const
{
  return (this == &right);
}

// --------------------------------------------------------------------
G4bool G4VProcess::operator!=(const G4VProcess& right) const
{
  return (this !=  &right);
}

// --------------------------------------------------------------------
void G4VProcess::DumpInfo() const
{
  G4cout << "Process Name " << theProcessName ;
  G4cout << " : Type[" << GetProcessTypeName(theProcessType) << "]";
  G4cout << " : SubType[" << theProcessSubType << "]"<< G4endl;
}

// --------------------------------------------------------------------
void G4VProcess::ProcessDescription(std::ostream& outFile) const
{
  outFile << "This process has not yet been described\n";
}
 
// --------------------------------------------------------------------
const G4String& G4VProcess::
GetPhysicsTableFileName( const G4ParticleDefinition* particle,
                         const G4String& directory,
                         const G4String& tableName,
                               G4bool ascii )
{
  G4String thePhysicsTableFileExt;
  if (ascii) thePhysicsTableFileExt = ".asc";
  else       thePhysicsTableFileExt = ".dat";

  thePhysicsTableFileName = directory + "/";
  thePhysicsTableFileName += tableName + "." +  theProcessName + ".";
  thePhysicsTableFileName += particle->GetParticleName()
                           + thePhysicsTableFileExt;
  
  return thePhysicsTableFileName;
}

// --------------------------------------------------------------------
void G4VProcess::BuildWorkerPhysicsTable(const G4ParticleDefinition& part)
{
  BuildPhysicsTable(part);
}

// --------------------------------------------------------------------
void G4VProcess::PrepareWorkerPhysicsTable(const G4ParticleDefinition& part)
{
  PreparePhysicsTable(part);
}

// --------------------------------------------------------------------
void G4VProcess::SetMasterProcess(G4VProcess* masterP)
{
  masterProcessShadow = masterP;
}
