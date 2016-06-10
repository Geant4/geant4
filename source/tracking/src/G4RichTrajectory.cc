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
// $Id: G4RichTrajectory.cc 91269 2015-06-29 07:05:59Z gcosmo $
//
// ---------------------------------------------------------------
//
// G4RichTrajectory.cc
//
// Contact:
//   Questions and comments on G4Trajectory, on which this is based,
//   should be sent to
//     Katsuya Amako  (e-mail: Katsuya.Amako@kek.jp)
//     Makoto  Asai   (e-mail: asai@kekvax.kek.jp)
//     Takashi Sasaki (e-mail: Takashi.Sasaki@kek.jp)
//   and on the extended code to:
//     John Allison   (e-mail: John.Allison@manchester.ac.uk)
//     Joseph Perl    (e-mail: perl@slac.stanford.edu)
//
// ---------------------------------------------------------------

#include "G4RichTrajectory.hh"
#include "G4RichTrajectoryPoint.hh"
#include "G4AttDefStore.hh"
#include "G4AttDef.hh"
#include "G4AttValue.hh"
#include "G4UIcommand.hh"
#include "G4UnitsTable.hh"
#include "G4VProcess.hh"

//#define G4ATTDEBUG
#ifdef G4ATTDEBUG
#include "G4AttCheck.hh"
#endif

#include <sstream>

G4ThreadLocal G4Allocator<G4RichTrajectory> *aRichTrajectoryAllocator = 0;

G4RichTrajectory::G4RichTrajectory():
  fpRichPointsContainer(0),
  fpCreatorProcess(0),
  fCreatorModelID(0),
  fpEndingProcess(0),
  fFinalKineticEnergy(0.)
{
}

G4RichTrajectory::G4RichTrajectory(const G4Track* aTrack):
  G4Trajectory(aTrack)  // Note: this initialises the base class data
			// members and, unfortunately but never mind,
			// creates a G4TrajectoryPoint in
			// TrajectoryPointContainer that we cannot
			// access because it's private.  We store the
			// same information (plus more) in a
			// G4RichTrajectoryPoint in the
			// RichTrajectoryPointsContainer
{
  fpInitialVolume = aTrack->GetTouchableHandle();
  fpInitialNextVolume = aTrack->GetNextTouchableHandle();
  fpCreatorProcess = aTrack->GetCreatorProcess();
  fCreatorModelID = aTrack->GetCreatorModelID();
  // On construction, set final values to initial values.
  // Final values are updated at the addition of every step - see AppendStep.
  fpFinalVolume = aTrack->GetTouchableHandle();
  fpFinalNextVolume = aTrack->GetNextTouchableHandle();
  fpEndingProcess = aTrack->GetCreatorProcess();
  fFinalKineticEnergy = aTrack->GetKineticEnergy();
  // Insert the first rich trajectory point (see note above)...
  fpRichPointsContainer = new RichTrajectoryPointsContainer;
  fpRichPointsContainer->push_back(new G4RichTrajectoryPoint(aTrack));
}

G4RichTrajectory::G4RichTrajectory(G4RichTrajectory & right):
  G4Trajectory(right)
{
  fpInitialVolume = right.fpInitialVolume;
  fpInitialNextVolume = right.fpInitialNextVolume;
  fpCreatorProcess = right.fpCreatorProcess;
  fCreatorModelID = right.fCreatorModelID;
  fpFinalVolume = right.fpFinalVolume;
  fpFinalNextVolume = right.fpFinalNextVolume;
  fpEndingProcess = right.fpEndingProcess;
  fFinalKineticEnergy = right.fFinalKineticEnergy;
  fpRichPointsContainer = new RichTrajectoryPointsContainer;
  for(size_t i=0;i<right.fpRichPointsContainer->size();i++)
  {
    G4RichTrajectoryPoint* rightPoint =
      (G4RichTrajectoryPoint*)((*(right.fpRichPointsContainer))[i]);
    fpRichPointsContainer->push_back(new G4RichTrajectoryPoint(*rightPoint));
  }
}

G4RichTrajectory::~G4RichTrajectory()
{
  if (fpRichPointsContainer) {
    //  fpRichPointsContainer->clearAndDestroy();
    size_t i;
    for(i=0;i<fpRichPointsContainer->size();i++){
      delete  (*fpRichPointsContainer)[i];
    }
    fpRichPointsContainer->clear();
    delete fpRichPointsContainer;
  }
}

void G4RichTrajectory::AppendStep(const G4Step* aStep)
{
  fpRichPointsContainer->push_back(new G4RichTrajectoryPoint(aStep));
  // Except for first step, which is a sort of virtual step to start
  // the track, compute the final values...
  const G4Track* track = aStep->GetTrack();
  const G4StepPoint* postStepPoint = aStep->GetPostStepPoint();
  if (track->GetCurrentStepNumber() > 0) {
    fpFinalVolume = track->GetTouchableHandle();
    fpFinalNextVolume = track->GetNextTouchableHandle();
    fpEndingProcess = postStepPoint->GetProcessDefinedStep();
    fFinalKineticEnergy =
      aStep->GetPreStepPoint()->GetKineticEnergy() -
      aStep->GetTotalEnergyDeposit();
  }
}
  
void G4RichTrajectory::MergeTrajectory(G4VTrajectory* secondTrajectory)
{
  if(!secondTrajectory) return;

  G4RichTrajectory* seco = (G4RichTrajectory*)secondTrajectory;
  G4int ent = seco->GetPointEntries();
  for(G4int i=1;i<ent;i++) {
    // initial point of the second trajectory should not be merged
    fpRichPointsContainer->push_back((*(seco->fpRichPointsContainer))[i]);
    //    fpRichPointsContainer->push_back(seco->fpRichPointsContainer->removeAt(1));
  }
  delete (*seco->fpRichPointsContainer)[0];
  seco->fpRichPointsContainer->clear();
}

void G4RichTrajectory::ShowTrajectory(std::ostream& os) const
{
  // Invoke the default implementation in G4VTrajectory...
  G4VTrajectory::ShowTrajectory(os);
  // ... or override with your own code here.
}

void G4RichTrajectory::DrawTrajectory() const
{
  // Invoke the default implementation in G4VTrajectory...
  G4VTrajectory::DrawTrajectory();
  // ... or override with your own code here.
}

const std::map<G4String,G4AttDef>* G4RichTrajectory::GetAttDefs() const
{
  G4bool isNew;
  std::map<G4String,G4AttDef>* store
    = G4AttDefStore::GetInstance("G4RichTrajectory",isNew);
  if (isNew) {

    // Get att defs from base class...
    *store = *(G4Trajectory::GetAttDefs());

    G4String ID;

    ID = "IVPath";
    (*store)[ID] = G4AttDef(ID,"Initial Volume Path",
			    "Physics","","G4String");

    ID = "INVPath";
    (*store)[ID] = G4AttDef(ID,"Initial Next Volume Path",
			    "Physics","","G4String");

    ID = "CPN";
    (*store)[ID] = G4AttDef(ID,"Creator Process Name",
                            "Physics","","G4String");

    ID = "CPTN";
    (*store)[ID] = G4AttDef(ID,"Creator Process Type Name",
                            "Physics","","G4String");
    
    ID = "CMID";
    (*store)[ID] = G4AttDef(ID,"Creator Model ID",
                            "Physics","","G4int");

    ID = "CMN";
    (*store)[ID] = G4AttDef(ID,"Creator Model Name",
                            "Physics","","G4String");

    ID = "FVPath";
    (*store)[ID] = G4AttDef(ID,"Final Volume Path",
			    "Physics","","G4String");

    ID = "FNVPath";
    (*store)[ID] = G4AttDef(ID,"Final Next Volume Path",
			    "Physics","","G4String");

    ID = "EPN";
    (*store)[ID] = G4AttDef(ID,"Ending Process Name",
			    "Physics","","G4String");

    ID = "EPTN";
    (*store)[ID] = G4AttDef(ID,"Ending Process Type Name",
			    "Physics","","G4String");

    ID = "FKE";
    (*store)[ID] = G4AttDef(ID,"Final kinetic energy",
			    "Physics","G4BestUnit","G4double");

  }

  return store;
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

std::vector<G4AttValue>* G4RichTrajectory::CreateAttValues() const
{
  // Create base class att values...
  std::vector<G4AttValue>* values = G4Trajectory::CreateAttValues();

  if (fpInitialVolume && fpInitialVolume->GetVolume()) {
    values->push_back(G4AttValue("IVPath",Path(fpInitialVolume),""));
  } else {
    values->push_back(G4AttValue("IVPath","None",""));
  }

  if (fpInitialNextVolume && fpInitialNextVolume->GetVolume()) {
    values->push_back(G4AttValue("INVPath",Path(fpInitialNextVolume),""));
  } else {
    values->push_back(G4AttValue("INVPath","None",""));
  }

  if (fpCreatorProcess) {
    values->push_back
    (G4AttValue("CPN",fpCreatorProcess->GetProcessName(),""));
    G4ProcessType type = fpCreatorProcess->GetProcessType();
    values->push_back
    (G4AttValue("CPTN",G4VProcess::GetProcessTypeName(type),""));
    values->push_back
    (G4AttValue("CMID",G4UIcommand::ConvertToString(fCreatorModelID),""));
    const G4String& creatorModelName =
    G4PhysicsModelCatalog::GetModelName(fCreatorModelID);
    values->push_back(G4AttValue("CMN",creatorModelName,""));
  } else {
    values->push_back(G4AttValue("CPN","None",""));
    values->push_back(G4AttValue("CPTN","None",""));
    values->push_back(G4AttValue("CMID","None",""));
    values->push_back(G4AttValue("CMN","None",""));
  }

  if (fpFinalVolume && fpFinalVolume->GetVolume()) {
    values->push_back(G4AttValue("FVPath",Path(fpFinalVolume),""));
  } else {
    values->push_back(G4AttValue("FVPath","None",""));
  }

  if (fpFinalNextVolume && fpFinalNextVolume->GetVolume()) {
    values->push_back(G4AttValue("FNVPath",Path(fpFinalNextVolume),""));
  } else {
    values->push_back(G4AttValue("FNVPath","None",""));
  }

  if (fpEndingProcess) {
    values->push_back(G4AttValue("EPN",fpEndingProcess->GetProcessName(),""));
    G4ProcessType type = fpEndingProcess->GetProcessType();
    values->push_back(G4AttValue("EPTN",G4VProcess::GetProcessTypeName(type),""));
  } else {
    values->push_back(G4AttValue("EPN","None",""));
    values->push_back(G4AttValue("EPTN","None",""));
  }

  values->push_back
    (G4AttValue("FKE",G4BestUnit(fFinalKineticEnergy,"Energy"),""));

#ifdef G4ATTDEBUG
  G4cout << G4AttCheck(values,GetAttDefs());
#endif

  return values;
}
