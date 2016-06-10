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
// $Id: G4RichTrajectoryPoint.cc 69003 2013-04-15 09:25:23Z gcosmo $
//
//
// ---------------------------------------------------------------
//
// G4RichTrajectoryPoint.cc
//
// Contact:
//   Questions and comments on G4TrajectoryPoint, on which this is based,
//   should be sent to
//     Katsuya Amako  (e-mail: Katsuya.Amako@kek.jp)
//     Makoto  Asai   (e-mail: asai@kekvax.kek.jp)
//     Takashi Sasaki (e-mail: Takashi.Sasaki@kek.jp)
//   and on the extended code to:
//     John Allison   (e-mail: John.Allison@manchester.ac.uk)
//     Joseph Perl    (e-mail: perl@slac.stanford.edu)
//
// ---------------------------------------------------------------

#include "G4RichTrajectoryPoint.hh"

#include "G4Track.hh"
#include "G4Step.hh"
#include "G4VProcess.hh"

#include "G4AttDefStore.hh"
#include "G4AttDef.hh"
#include "G4AttValue.hh"
#include "G4UnitsTable.hh"

//#define G4ATTDEBUG
#ifdef G4ATTDEBUG
#include "G4AttCheck.hh"
#endif

#include <sstream>

G4ThreadLocal G4Allocator<G4RichTrajectoryPoint> *aRichTrajectoryPointAllocator = 0;

G4RichTrajectoryPoint::G4RichTrajectoryPoint():
  fpAuxiliaryPointVector(0),
  fTotEDep(0.),
  fRemainingEnergy(0.),
  fpProcess(0),
  fPreStepPointStatus(fUndefined),
  fPostStepPointStatus(fUndefined),
  fPreStepPointGlobalTime(0),
  fPostStepPointGlobalTime(0),
  fPreStepPointWeight(1.),
  fPostStepPointWeight(1.)
{}

G4RichTrajectoryPoint::G4RichTrajectoryPoint(const G4Track* aTrack):
  G4TrajectoryPoint(aTrack->GetPosition()),
  fpAuxiliaryPointVector(0),
  fTotEDep(0.),
  fRemainingEnergy(aTrack->GetKineticEnergy()),
  fpProcess(0),
  fPreStepPointStatus(fUndefined),
  fPostStepPointStatus(fUndefined),
  fPreStepPointGlobalTime(aTrack->GetGlobalTime()),
  fPostStepPointGlobalTime(aTrack->GetGlobalTime()),
  fpPreStepPointVolume(aTrack->GetTouchableHandle()),
  fpPostStepPointVolume(aTrack->GetNextTouchableHandle()),
  fPreStepPointWeight(aTrack->GetWeight()),
  fPostStepPointWeight(aTrack->GetWeight())
{}

G4RichTrajectoryPoint::G4RichTrajectoryPoint(const G4Step* aStep):
  G4TrajectoryPoint(aStep->GetPostStepPoint()->GetPosition()),
  fpAuxiliaryPointVector(aStep->GetPointerToVectorOfAuxiliaryPoints()),
  fTotEDep(aStep->GetTotalEnergyDeposit())
{
  G4StepPoint* preStepPoint = aStep->GetPreStepPoint();
  G4StepPoint* postStepPoint = aStep->GetPostStepPoint();
  if (aStep->GetTrack()->GetCurrentStepNumber() <= 0) {  // First step
    fRemainingEnergy = aStep->GetTrack()->GetKineticEnergy();
  } else {
    fRemainingEnergy = preStepPoint->GetKineticEnergy() - fTotEDep;
  }
  fpProcess = postStepPoint->GetProcessDefinedStep();
  fPreStepPointStatus = preStepPoint->GetStepStatus();
  fPostStepPointStatus = postStepPoint->GetStepStatus();
  fPreStepPointGlobalTime = preStepPoint->GetGlobalTime();
  fPostStepPointGlobalTime = postStepPoint->GetGlobalTime();
  fpPreStepPointVolume = preStepPoint->GetTouchableHandle();
  fpPostStepPointVolume = postStepPoint->GetTouchableHandle();
  fPreStepPointWeight = preStepPoint->GetWeight();
  fPostStepPointWeight = postStepPoint->GetWeight();
}

G4RichTrajectoryPoint::G4RichTrajectoryPoint
(const G4RichTrajectoryPoint &right):
  G4TrajectoryPoint(right),
  fpAuxiliaryPointVector(right.fpAuxiliaryPointVector),
  fTotEDep(right.fTotEDep),
  fRemainingEnergy(right.fRemainingEnergy),
  fpProcess(right.fpProcess),
  fPreStepPointStatus(right.fPreStepPointStatus),
  fPostStepPointStatus(right.fPostStepPointStatus),
  fPreStepPointGlobalTime(right.fPreStepPointGlobalTime),
  fPostStepPointGlobalTime(right.fPostStepPointGlobalTime),
  fpPreStepPointVolume(right.fpPreStepPointVolume),
  fpPostStepPointVolume(right.fpPostStepPointVolume),
  fPreStepPointWeight(right.fPreStepPointWeight),
  fPostStepPointWeight(right.fPostStepPointWeight)
{}

G4RichTrajectoryPoint::~G4RichTrajectoryPoint()
{
  if(fpAuxiliaryPointVector) {
    /*
    G4cout << "Deleting fpAuxiliaryPointVector at "
	   << (void*) fpAuxiliaryPointVector
	   << G4endl;
    */
    delete fpAuxiliaryPointVector;
  }
}

const std::map<G4String,G4AttDef>*
G4RichTrajectoryPoint::GetAttDefs() const
{
  G4bool isNew;
  std::map<G4String,G4AttDef>* store
    = G4AttDefStore::GetInstance("G4RichTrajectoryPoint",isNew);
  if (isNew) {

    // Copy base class att defs...
    *store = *(G4TrajectoryPoint::GetAttDefs());

    G4String ID;

    ID = "Aux";
    (*store)[ID] = G4AttDef(ID, "Auxiliary Point Position",
			    "Physics","G4BestUnit","G4ThreeVector");
    ID = "TED";
    (*store)[ID] = G4AttDef(ID,"Total Energy Deposit",
			    "Physics","G4BestUnit","G4double");
    ID = "RE";
    (*store)[ID] = G4AttDef(ID,"Remaining Energy",
			    "Physics","G4BestUnit","G4double");
    ID = "PDS";
    (*store)[ID] = G4AttDef(ID,"Process Defined Step",
			    "Physics","","G4String");
    ID = "PTDS";
    (*store)[ID] = G4AttDef(ID,"Process Type Defined Step",
			    "Physics","","G4String");
    ID = "PreStatus";
    (*store)[ID] = G4AttDef(ID,"Pre-step-point status",
			    "Physics","","G4String");
    ID = "PostStatus";
    (*store)[ID] = G4AttDef(ID,"Post-step-point status",
			    "Physics","","G4String");
    ID = "PreT";
    (*store)[ID] = G4AttDef(ID,"Pre-step-point global time",
			    "Physics","G4BestUnit","G4double");
    ID = "PostT";
    (*store)[ID] = G4AttDef(ID,"Post-step-point global time",
			    "Physics","G4BestUnit","G4double");
    ID = "PreVPath";
    (*store)[ID] = G4AttDef(ID,"Pre-step Volume Path",
                            "Physics","","G4String");
    ID = "PostVPath";
    (*store)[ID] = G4AttDef(ID,"Post-step Volume Path",
                            "Physics","","G4String");
    ID = "PreW";
    (*store)[ID] = G4AttDef(ID,"Pre-step-point weight",
			    "Physics","","G4double");
    ID = "PostW";
    (*store)[ID] = G4AttDef(ID,"Post-step-point weight",
			    "Physics","","G4double");
  }
  return store;
}

static G4String Status(G4StepStatus stps)
{
  G4String status;
  switch (stps) {
  case fWorldBoundary:         status = "fWorldBoundary"; break;
  case fGeomBoundary:          status = "fGeomBoundary"; break;
  case fAtRestDoItProc:        status = "fAtRestDoItProc"; break;
  case fAlongStepDoItProc:     status = "fAlongStepDoItProc"; break;
  case fPostStepDoItProc:      status = "fPostStepDoItProc"; break;
  case fUserDefinedLimit:      status = "fUserDefinedLimit"; break;
  case fExclusivelyForcedProc: status = "fExclusivelyForcedProc"; break;
  case fUndefined:             status = "fUndefined"; break;
  default:                     status = "Not recognised"; break;
  }
  return status;
}

static G4String Path(const G4TouchableHandle& th)
{
  std::ostringstream oss;
  G4int depth = th->GetHistoryDepth();
  for (G4int i = depth; i >= 0; --i) {
    oss << th->GetVolume(i)->GetName()
	<< ':' << th->GetCopyNumber(i);
    if (i != 0) oss << '/';
  }
  return oss.str();
}

std::vector<G4AttValue>* G4RichTrajectoryPoint::CreateAttValues() const
{
  // Create base class att values...
  std::vector<G4AttValue>* values = G4TrajectoryPoint::CreateAttValues();

  if (fpAuxiliaryPointVector) {
    std::vector<G4ThreeVector>::iterator iAux;
    for (iAux = fpAuxiliaryPointVector->begin();
	 iAux != fpAuxiliaryPointVector->end(); ++iAux) {
      values->push_back(G4AttValue("Aux",G4BestUnit(*iAux,"Length"),""));
    }
  }

  values->push_back(G4AttValue("TED",G4BestUnit(fTotEDep,"Energy"),""));

  values->push_back(G4AttValue("RE",G4BestUnit(fRemainingEnergy,"Energy"),""));

  if (fpProcess) {
    values->push_back
      (G4AttValue("PDS",fpProcess->GetProcessName(),""));
    values->push_back
      (G4AttValue
       ("PTDS",G4VProcess::GetProcessTypeName(fpProcess->GetProcessType()),
	""));
  } else {
    values->push_back(G4AttValue("PDS","None",""));
    values->push_back(G4AttValue("PTDS","None",""));
  }

  values->push_back
    (G4AttValue("PreStatus",Status(fPreStepPointStatus),""));

  values->push_back
    (G4AttValue("PostStatus",Status(fPostStepPointStatus),""));

  values->push_back
    (G4AttValue("PreT",G4BestUnit(fPreStepPointGlobalTime,"Time"),""));

  values->push_back
    (G4AttValue("PostT",G4BestUnit(fPostStepPointGlobalTime,"Time"),""));

  if (fpPreStepPointVolume && fpPreStepPointVolume->GetVolume()) {
    values->push_back(G4AttValue("PreVPath",Path(fpPreStepPointVolume),""));
  } else {
    values->push_back(G4AttValue("PreVPath","None",""));
  }

  if (fpPostStepPointVolume && fpPostStepPointVolume->GetVolume()) {
    values->push_back(G4AttValue("PostVPath",Path(fpPostStepPointVolume),""));
  } else {
    values->push_back(G4AttValue("PostVPath","None",""));
  }

  {
    std::ostringstream oss;
    oss << fPreStepPointWeight;
    values->push_back
      (G4AttValue("PreW",oss.str(),""));
  }

  {
    std::ostringstream oss;
    oss << fPostStepPointWeight;
    values->push_back
      (G4AttValue("PostW",oss.str(),""));
  }

#ifdef G4ATTDEBUG
  G4cout << G4AttCheck(values,GetAttDefs());
#endif

  return values;
}
