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
// $Id: G4RichTrajectory.cc,v 1.6 2006/10/16 13:43:43 allison Exp $
// GEANT4 tag $Name: geant4-09-01 $
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
#include "G4VProcess.hh"

//#define G4ATTDEBUG
#ifdef G4ATTDEBUG
#include "G4AttCheck.hh"
#endif

G4Allocator<G4RichTrajectory> aRichTrajectoryAllocator;

G4RichTrajectory::G4RichTrajectory():
  fpRichPointsContainer(0),
  fpInitialVolume(0),
  fpInitialNextVolume(0),
  fpCreatorProcess(0)
{}

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
  fpInitialVolume = aTrack->GetVolume();
  fpInitialNextVolume = aTrack->GetNextVolume();
  fpCreatorProcess = aTrack->GetCreatorProcess();
  fpRichPointsContainer = new RichTrajectoryPointsContainer;
  // Insert the first rich trajectory point (see note above)...
  fpRichPointsContainer->push_back(new G4RichTrajectoryPoint(aTrack));
}

G4RichTrajectory::G4RichTrajectory(G4RichTrajectory & right):
  G4Trajectory(right)
{
  fpInitialVolume = right.fpInitialVolume;
  fpInitialNextVolume = right.fpInitialNextVolume;
  fpCreatorProcess = right.fpCreatorProcess;
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
}
  
void G4RichTrajectory::MergeTrajectory(G4VTrajectory* secondTrajectory)
{
  if(!secondTrajectory) return;

  G4Trajectory::MergeTrajectory(secondTrajectory);

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

const std::map<G4String,G4AttDef>* G4RichTrajectory::GetAttDefs() const
{
  G4bool isNew;
  std::map<G4String,G4AttDef>* store
    = G4AttDefStore::GetInstance("G4RichTrajectory",isNew);
  if (isNew) {

    // Get att defs from base class...
    *store = *(G4Trajectory::GetAttDefs());

    G4String ID;

    ID = "IVN";
    (*store)[ID] = G4AttDef(ID,"Initial Volume Name",
			    "Physics","","G4String");

    ID = "INVN";
    (*store)[ID] = G4AttDef(ID,"Initial Next Volume Name",
			    "Physics","","G4String");

    ID = "CPN";
    (*store)[ID] = G4AttDef(ID,"Creator Process Name",
			    "Physics","","G4String");

    ID = "CPTN";
    (*store)[ID] = G4AttDef(ID,"Creator Process Type Name",
			    "Physics","","G4String");

  }

  return store;
}

std::vector<G4AttValue>* G4RichTrajectory::CreateAttValues() const
{
  // Create base class att values...
  std::vector<G4AttValue>* values = G4Trajectory::CreateAttValues();

  values->push_back(G4AttValue("IVN",fpInitialVolume->GetName(),""));

  values->push_back(G4AttValue("INVN",fpInitialNextVolume->GetName(),""));

  if (fpCreatorProcess) {
    values->push_back(G4AttValue("CPN",fpCreatorProcess->GetProcessName(),""));
    G4ProcessType type = fpCreatorProcess->GetProcessType();
    values->push_back(G4AttValue("CPTN",G4VProcess::GetProcessTypeName(type),""));
  } else {
    values->push_back(G4AttValue("CPN","User Defined",""));
    values->push_back(G4AttValue("CPTN","User",""));
  }

#ifdef G4ATTDEBUG
  G4cout << G4AttCheck(values,GetAttDefs());
#endif

  return values;
}
