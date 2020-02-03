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
/// \file runAndEvent/RE04/src/RE04Trajectory.cc
/// \brief Implementation of the RE04Trajectory class
//
//
#include "RE04Trajectory.hh"
#include "RE04TrajectoryPoint.hh"
#include "G4ParticleTable.hh"
#include "G4AttDefStore.hh"
#include "G4AttDef.hh"
#include "G4AttValue.hh"
#include "G4UIcommand.hh"
#include "G4UnitsTable.hh"

//#define G4ATTDEBUG
#ifdef G4ATTDEBUG
#include "G4AttCheck.hh"
#endif

G4ThreadLocal G4Allocator<RE04Trajectory> * faTrajAllocator = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RE04Trajectory::RE04Trajectory()
:  G4VTrajectory(),
   fPositionRecord(0), fTrackID(0), fParentID(0),
   fPDGEncoding( 0 ), fPDGCharge(0.0), fParticleName(""),
   fInitialKineticEnergy( 0. ), fInitialMomentum( G4ThreeVector() )
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RE04Trajectory::RE04Trajectory(const G4Track* aTrack)
{
   G4ParticleDefinition * fpParticleDefinition = aTrack->GetDefinition();
   fParticleName = fpParticleDefinition->GetParticleName();
   fPDGCharge = fpParticleDefinition->GetPDGCharge();
   fPDGEncoding = fpParticleDefinition->GetPDGEncoding();
   fTrackID = aTrack->GetTrackID();
   fParentID = aTrack->GetParentID();
   fInitialKineticEnergy = aTrack->GetKineticEnergy();
   fInitialMomentum = aTrack->GetMomentum();
   fPositionRecord = new TrajectoryPointContainer();
   // Following is for the first trajectory point
   fPositionRecord->push_back(new RE04TrajectoryPoint(
      aTrack->GetPosition(),aTrack->GetMaterial()));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RE04Trajectory::RE04Trajectory(RE04Trajectory & right):G4VTrajectory()
{
  fParticleName = right.fParticleName;
  fPDGCharge = right.fPDGCharge;
  fPDGEncoding = right.fPDGEncoding;
  fTrackID = right.fTrackID;
  fParentID = right.fParentID;
  fInitialKineticEnergy = right.fInitialKineticEnergy;
  fInitialMomentum = right.fInitialMomentum;
  fPositionRecord = new TrajectoryPointContainer();

  for(size_t i=0;i<right.fPositionRecord->size();i++)
  {
    RE04TrajectoryPoint* rightPoint
      = (RE04TrajectoryPoint*)((*(right.fPositionRecord))[i]);
    fPositionRecord->push_back(new RE04TrajectoryPoint(*rightPoint));
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RE04Trajectory::~RE04Trajectory()
{
  if (fPositionRecord) {
    //  fPositionRecord->clearAndDestroy();
    size_t i;
    for(i=0;i<fPositionRecord->size();i++){
      delete  (*fPositionRecord)[i];
    }
    fPositionRecord->clear();
    delete fPositionRecord;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RE04Trajectory::ShowTrajectory(std::ostream& os) const
{
  // Invoke the default implementation in G4VTrajectory...
  G4VTrajectory::ShowTrajectory(os);
  // ... or override with your own code here.
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RE04Trajectory::DrawTrajectory() const
{
  // Invoke the default implementation in G4VTrajectory...
  G4VTrajectory::DrawTrajectory();
  // ... or override with your own code here.
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
const std::map<G4String,G4AttDef>* RE04Trajectory::GetAttDefs() const
{
  G4bool isNew;
  std::map<G4String,G4AttDef>* store
    = G4AttDefStore::GetInstance("RE04Trajectory",isNew);
  if (isNew) {

    G4String id("ID");
    (*store)[id] = G4AttDef(id,"Track ID","Physics","","G4int");

    G4String pid("PID");
    (*store)[pid] = G4AttDef(pid,"Parent ID","Physics","","G4int");

    G4String pn("PN");
    (*store)[pn] = G4AttDef(pn,"Particle Name","Physics","","G4String");

    G4String ch("Ch");
    (*store)[ch] = G4AttDef(ch,"Charge","Physics","e+","G4double");

    G4String pdg("PDG");
    (*store)[pdg] = G4AttDef(pdg,"PDG Encoding","Physics","","G4int");

    G4String ike("IKE");
    (*store)[ike] = 
      G4AttDef(ike, "Initial kinetic energy",
               "Physics","G4BestUnit","G4double");

    G4String iMom("IMom");
    (*store)[iMom] = G4AttDef(iMom, "Initial momentum",
                              "Physics","G4BestUnit","G4ThreeVector");

    G4String iMag("IMag");
    (*store)[iMag] = 
      G4AttDef(iMag, "Initial momentum magnitude",
               "Physics","G4BestUnit","G4double");

    G4String ntp("NTP");
    (*store)[ntp] = G4AttDef(ntp,"No. of points","Physics","","G4int");

  }
  return store;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
std::vector<G4AttValue>* RE04Trajectory::CreateAttValues() const
{
  std::vector<G4AttValue>* values = new std::vector<G4AttValue>;

  values->push_back
    (G4AttValue("ID",G4UIcommand::ConvertToString(fTrackID),""));

  values->push_back
    (G4AttValue("PID",G4UIcommand::ConvertToString(fParentID),""));

  values->push_back(G4AttValue("PN",fParticleName,""));

  values->push_back
    (G4AttValue("Ch",G4UIcommand::ConvertToString(fPDGCharge),""));

  values->push_back
    (G4AttValue("PDG",G4UIcommand::ConvertToString(fPDGEncoding),""));

  values->push_back
    (G4AttValue("IKE",G4BestUnit(fInitialKineticEnergy,"Energy"),""));

  values->push_back
    (G4AttValue("IMom",G4BestUnit(fInitialMomentum,"Energy"),""));

  values->push_back
    (G4AttValue("IMag",G4BestUnit(fInitialMomentum.mag(),"Energy"),""));

  values->push_back
    (G4AttValue("NTP",G4UIcommand::ConvertToString(GetPointEntries()),""));

#ifdef G4ATTDEBUG
  G4cout << G4AttCheck(values,GetAttDefs());
#endif

  return values;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RE04Trajectory::AppendStep(const G4Step* aStep)
{
   fPositionRecord->push_back( new RE04TrajectoryPoint(
        aStep->GetPostStepPoint()->GetPosition(),
        aStep->GetPreStepPoint()->GetMaterial() ));
}
  
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4ParticleDefinition* RE04Trajectory::GetParticleDefinition()
{
   return (G4ParticleTable::GetParticleTable()->FindParticle(fParticleName));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RE04Trajectory::MergeTrajectory(G4VTrajectory* secondTrajectory)
{
  if(!secondTrajectory) return;

  RE04Trajectory* seco = (RE04Trajectory*)secondTrajectory;
  G4int ent = seco->GetPointEntries();
  for(G4int i=1;i<ent;i++) // initial point of the second trajectory
                           // should not be merged
  { 
    fPositionRecord->push_back((*(seco->fPositionRecord))[i]);
    //    fPositionRecord->push_back(seco->fPositionRecord->removeAt(1));
  }
  delete (*seco->fPositionRecord)[0];
  seco->fPositionRecord->clear();
}


