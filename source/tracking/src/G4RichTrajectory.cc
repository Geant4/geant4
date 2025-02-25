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
// G4RichTrajectory class implementation
//
// Contact:
//   Questions and comments on G4Trajectory, on which this is based,
//   should be sent to
//     Katsuya Amako  (e-mail: Katsuya.Amako@kek.jp)
//     Makoto  Asai   (e-mail: asai@slac.stanford.edu)
//     Takashi Sasaki (e-mail: Takashi.Sasaki@kek.jp)
//   and on the extended code to:
//     John Allison   (e-mail: John.Allison@manchester.ac.uk)
//     Joseph Perl    (e-mail: perl@slac.stanford.edu)
// --------------------------------------------------------------------

#include "G4RichTrajectory.hh"
#include "G4ClonedRichTrajectory.hh"

#include "G4ParticleTable.hh"
#include "G4AttDef.hh"
#include "G4AttDefStore.hh"
#include "G4AttValue.hh"
#include "G4PhysicsModelCatalog.hh"
#include "G4RichTrajectoryPoint.hh"
#include "G4UIcommand.hh"
#include "G4UnitsTable.hh"
#include "G4VProcess.hh"

namespace {
 G4Mutex CloneRichTrajectoryMutex = G4MUTEX_INITIALIZER;
}

// #define G4ATTDEBUG
#ifdef G4ATTDEBUG
#  include "G4AttCheck.hh"
#endif

#include <sstream>

G4Allocator<G4RichTrajectory>*& aRichTrajectoryAllocator()
{
  G4ThreadLocalStatic G4Allocator<G4RichTrajectory>* _instance = nullptr;
  return _instance;
}

G4RichTrajectory::G4RichTrajectory(const G4Track* aTrack)
{
  G4ParticleDefinition* fpParticleDefinition = aTrack->GetDefinition();
  ParticleName = fpParticleDefinition->GetParticleName();
  PDGCharge = fpParticleDefinition->GetPDGCharge();
  PDGEncoding = fpParticleDefinition->GetPDGEncoding();
  fTrackID = aTrack->GetTrackID();
  fParentID = aTrack->GetParentID();
  initialKineticEnergy = aTrack->GetKineticEnergy();
  initialMomentum = aTrack->GetMomentum();
  positionRecord = new G4TrajectoryPointContainer();

  // Following is for the first trajectory point
  positionRecord->push_back(new G4RichTrajectoryPoint(aTrack));

  fpInitialVolume = aTrack->GetTouchableHandle();
  fpInitialNextVolume = aTrack->GetNextTouchableHandle();
  fpCreatorProcess = aTrack->GetCreatorProcess();
  fCreatorModelID = aTrack->GetCreatorModelID();

  // On construction, set final values to initial values.
  // Final values are updated at the addition of every step - see AppendStep.
  //
  fpFinalVolume = aTrack->GetTouchableHandle();
  fpFinalNextVolume = aTrack->GetNextTouchableHandle();
  fpEndingProcess = aTrack->GetCreatorProcess();
  fFinalKineticEnergy = aTrack->GetKineticEnergy();

  // Insert the first rich trajectory point (see note above)...
  //
  fpRichPointContainer = new G4TrajectoryPointContainer;
  fpRichPointContainer->push_back(new G4RichTrajectoryPoint(aTrack));
}

G4RichTrajectory::G4RichTrajectory(G4RichTrajectory& right) 
{
  ParticleName = right.ParticleName;
  PDGCharge = right.PDGCharge;
  PDGEncoding = right.PDGEncoding;
  fTrackID = right.fTrackID;
  fParentID = right.fParentID;
  initialKineticEnergy = right.initialKineticEnergy;
  initialMomentum = right.initialMomentum;
  positionRecord = new G4TrajectoryPointContainer();

  for (auto& i : *right.positionRecord) {
    auto rightPoint = (G4RichTrajectoryPoint*)i;
    positionRecord->push_back(new G4RichTrajectoryPoint(*rightPoint));
  }

  fpInitialVolume = right.fpInitialVolume;
  fpInitialNextVolume = right.fpInitialNextVolume;
  fpCreatorProcess = right.fpCreatorProcess;
  fCreatorModelID = right.fCreatorModelID;
  fpFinalVolume = right.fpFinalVolume;
  fpFinalNextVolume = right.fpFinalNextVolume;
  fpEndingProcess = right.fpEndingProcess;
  fFinalKineticEnergy = right.fFinalKineticEnergy;
  fpRichPointContainer = new G4TrajectoryPointContainer;
  for (auto& i : *right.fpRichPointContainer) {
    auto rightPoint = (G4RichTrajectoryPoint*)i;
    fpRichPointContainer->push_back(new G4RichTrajectoryPoint(*rightPoint));
  }
}

G4RichTrajectory::~G4RichTrajectory()
{
  if (fpRichPointContainer != nullptr) {
    for (auto& i : *fpRichPointContainer) {
      delete i;
    }
    fpRichPointContainer->clear();
    delete fpRichPointContainer;
  }
}

void G4RichTrajectory::AppendStep(const G4Step* aStep)
{
  fpRichPointContainer->push_back(new G4RichTrajectoryPoint(aStep));

  // Except for first step, which is a sort of virtual step to start
  // the track, compute the final values...
  //
  const G4Track* track = aStep->GetTrack();
  const G4StepPoint* postStepPoint = aStep->GetPostStepPoint();
  if (track->GetCurrentStepNumber() > 0) {
    fpFinalVolume = track->GetTouchableHandle();
    fpFinalNextVolume = track->GetNextTouchableHandle();
    fpEndingProcess = postStepPoint->GetProcessDefinedStep();
    fFinalKineticEnergy =
      aStep->GetPreStepPoint()->GetKineticEnergy() - aStep->GetTotalEnergyDeposit();
  }
}

void G4RichTrajectory::MergeTrajectory(G4VTrajectory* secondTrajectory)
{
  if (secondTrajectory == nullptr) return;

  auto seco = (G4RichTrajectory*)secondTrajectory;
  G4int ent = seco->GetPointEntries();
  for (G4int i = 1; i < ent; ++i) {
    // initial point of the second trajectory should not be merged
    //
    fpRichPointContainer->push_back((*(seco->fpRichPointContainer))[i]);
  }
  delete (*seco->fpRichPointContainer)[0];
  seco->fpRichPointContainer->clear();
}

void G4RichTrajectory::ShowTrajectory(std::ostream& os) const
{
  // Invoke the default implementation in G4VTrajectory...
  //
  G4VTrajectory::ShowTrajectory(os);

  // ... or override with your own code here.
}

void G4RichTrajectory::DrawTrajectory() const
{
  // Invoke the default implementation in G4VTrajectory...
  //
  G4VTrajectory::DrawTrajectory();

  // ... or override with your own code here.
}

const std::map<G4String, G4AttDef>* G4RichTrajectory::GetAttDefs() const
{
  G4bool isNew;
  std::map<G4String, G4AttDef>* store = G4AttDefStore::GetInstance("G4RichTrajectory", isNew);
  if (isNew) {
    G4String ID;

    ID = "ID";
    (*store)[ID] = G4AttDef(ID, "Track ID", "Physics", "", "G4int");

    ID = "PID";
    (*store)[ID] = G4AttDef(ID, "Parent ID", "Physics", "", "G4int");

    ID = "PN";
    (*store)[ID] = G4AttDef(ID, "Particle Name", "Physics", "", "G4String");

    ID = "Ch";
    (*store)[ID] = G4AttDef(ID, "Charge", "Physics", "e+", "G4double");

    ID = "PDG";
    (*store)[ID] = G4AttDef(ID, "PDG Encoding", "Physics", "", "G4int");

    ID = "IKE";
    (*store)[ID] = G4AttDef(ID, "Initial kinetic energy", "Physics", "G4BestUnit", "G4double");

    ID = "IMom";
    (*store)[ID] = G4AttDef(ID, "Initial momentum", "Physics", "G4BestUnit", "G4ThreeVector");

    ID = "IMag";
    (*store)[ID] = G4AttDef(ID, "Initial momentum magnitude", "Physics", "G4BestUnit", "G4double");

    ID = "NTP";
    (*store)[ID] = G4AttDef(ID, "No. of points", "Physics", "", "G4int");

    ID = "IVPath";
    (*store)[ID] = G4AttDef(ID, "Initial Volume Path", "Physics", "", "G4String");

    ID = "INVPath";
    (*store)[ID] = G4AttDef(ID, "Initial Next Volume Path", "Physics", "", "G4String");

    ID = "CPN";
    (*store)[ID] = G4AttDef(ID, "Creator Process Name", "Physics", "", "G4String");

    ID = "CPTN";
    (*store)[ID] = G4AttDef(ID, "Creator Process Type Name", "Physics", "", "G4String");

    ID = "CMID";
    (*store)[ID] = G4AttDef(ID, "Creator Model ID", "Physics", "", "G4int");

    ID = "CMN";
    (*store)[ID] = G4AttDef(ID, "Creator Model Name", "Physics", "", "G4String");

    ID = "FVPath";
    (*store)[ID] = G4AttDef(ID, "Final Volume Path", "Physics", "", "G4String");

    ID = "FNVPath";
    (*store)[ID] = G4AttDef(ID, "Final Next Volume Path", "Physics", "", "G4String");

    ID = "EPN";
    (*store)[ID] = G4AttDef(ID, "Ending Process Name", "Physics", "", "G4String");

    ID = "EPTN";
    (*store)[ID] = G4AttDef(ID, "Ending Process Type Name", "Physics", "", "G4String");

    ID = "FKE";
    (*store)[ID] = G4AttDef(ID, "Final kinetic energy", "Physics", "G4BestUnit", "G4double");
  }

  return store;
}

static G4String Path(const G4TouchableHandle& th)
{
  std::ostringstream oss;
  G4int depth = th->GetHistoryDepth();
  for (G4int i = depth; i >= 0; --i) {
    oss << th->GetVolume(i)->GetName() << ':' << th->GetCopyNumber(i);
    if (i != 0) oss << '/';
  }
  return oss.str();
}

std::vector<G4AttValue>* G4RichTrajectory::CreateAttValues() const
{
  // Create base class att values...
  //std::vector<G4AttValue>* values = G4VTrajectory::CreateAttValues();
  auto values = new std::vector<G4AttValue>;
  values->push_back(G4AttValue("ID", G4UIcommand::ConvertToString(fTrackID), ""));
  values->push_back(G4AttValue("PID", G4UIcommand::ConvertToString(fParentID), ""));
  values->push_back(G4AttValue("PN", ParticleName, ""));
  values->push_back(G4AttValue("Ch", G4UIcommand::ConvertToString(PDGCharge), ""));
  values->push_back(G4AttValue("PDG", G4UIcommand::ConvertToString(PDGEncoding), ""));
  values->push_back(G4AttValue("IKE", G4BestUnit(initialKineticEnergy, "Energy"), ""));
  values->push_back(G4AttValue("IMom", G4BestUnit(initialMomentum, "Energy"), ""));
  values->push_back(G4AttValue("IMag", G4BestUnit(initialMomentum.mag(), "Energy"), ""));
  values->push_back(G4AttValue("NTP", G4UIcommand::ConvertToString(GetPointEntries()), ""));

  if (fpInitialVolume && (fpInitialVolume->GetVolume() != nullptr)) {
    values->push_back(G4AttValue("IVPath", Path(fpInitialVolume), ""));
  }
  else {
    values->push_back(G4AttValue("IVPath", "None", ""));
  }

  if (fpInitialNextVolume && (fpInitialNextVolume->GetVolume() != nullptr)) {
    values->push_back(G4AttValue("INVPath", Path(fpInitialNextVolume), ""));
  }
  else {
    values->push_back(G4AttValue("INVPath", "None", ""));
  }

  if (fpCreatorProcess != nullptr) {
    values->push_back(G4AttValue("CPN", fpCreatorProcess->GetProcessName(), ""));
    G4ProcessType type = fpCreatorProcess->GetProcessType();
    values->push_back(G4AttValue("CPTN", G4VProcess::GetProcessTypeName(type), ""));
    values->push_back(G4AttValue("CMID", G4UIcommand::ConvertToString(fCreatorModelID), ""));
    const G4String& creatorModelName = G4PhysicsModelCatalog::GetModelNameFromID(fCreatorModelID);
    values->push_back(G4AttValue("CMN", creatorModelName, ""));
  }
  else {
    values->push_back(G4AttValue("CPN", "None", ""));
    values->push_back(G4AttValue("CPTN", "None", ""));
    values->push_back(G4AttValue("CMID", "None", ""));
    values->push_back(G4AttValue("CMN", "None", ""));
  }

  if (fpFinalVolume && (fpFinalVolume->GetVolume() != nullptr)) {
    values->push_back(G4AttValue("FVPath", Path(fpFinalVolume), ""));
  }
  else {
    values->push_back(G4AttValue("FVPath", "None", ""));
  }

  if (fpFinalNextVolume && (fpFinalNextVolume->GetVolume() != nullptr)) {
    values->push_back(G4AttValue("FNVPath", Path(fpFinalNextVolume), ""));
  }
  else {
    values->push_back(G4AttValue("FNVPath", "None", ""));
  }

  if (fpEndingProcess != nullptr) {
    values->push_back(G4AttValue("EPN", fpEndingProcess->GetProcessName(), ""));
    G4ProcessType type = fpEndingProcess->GetProcessType();
    values->push_back(G4AttValue("EPTN", G4VProcess::GetProcessTypeName(type), ""));
  }
  else {
    values->push_back(G4AttValue("EPN", "None", ""));
    values->push_back(G4AttValue("EPTN", "None", ""));
  }

  values->push_back(G4AttValue("FKE", G4BestUnit(fFinalKineticEnergy, "Energy"), ""));

#ifdef G4ATTDEBUG
  G4cout << G4AttCheck(values, GetAttDefs());
#endif

  return values;
}

G4ParticleDefinition* G4RichTrajectory::GetParticleDefinition()
{
  return (G4ParticleTable::GetParticleTable()->FindParticle(ParticleName));
}

G4VTrajectory* G4RichTrajectory::CloneForMaster() const
{
  G4AutoLock lock(&CloneRichTrajectoryMutex);
  auto* cloned = new G4ClonedRichTrajectory(*this);
  return cloned;
}

