//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4RichTrajectoryPoint.cc,v 1.1 2005/11/24 12:47:36 allison Exp $
// GEANT4 tag $Name: geant4-08-00 $
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

G4Allocator<G4RichTrajectoryPoint> aRichTrajectoryPointAllocator;

G4RichTrajectoryPoint::G4RichTrajectoryPoint():
  fTotEDep(0.),
  fpProcess(0)
{}

G4RichTrajectoryPoint::G4RichTrajectoryPoint(const G4Track* aTrack):
  G4TrajectoryPoint(aTrack->GetPosition()),
  fTotEDep(0.),
  fpProcess(0)
{}

G4RichTrajectoryPoint::G4RichTrajectoryPoint(const G4Step* aStep):
  G4TrajectoryPoint(aStep->GetPostStepPoint()->GetPosition()),
  fTotEDep(aStep->GetTotalEnergyDeposit()),
  fpProcess(aStep->GetPostStepPoint()->GetProcessDefinedStep())
{}

G4RichTrajectoryPoint::G4RichTrajectoryPoint
(const G4RichTrajectoryPoint &right):
  G4TrajectoryPoint(right)
{}

G4RichTrajectoryPoint::~G4RichTrajectoryPoint()
{}

const std::map<G4String,G4AttDef>*
G4RichTrajectoryPoint::GetAttDefs() const
{
  G4bool isNew;
  std::map<G4String,G4AttDef>* store
    = G4AttDefStore::GetInstance("G4RichTrajectoryPoint",isNew);
  if (isNew) {
    *store = *(G4TrajectoryPoint::GetAttDefs());

    G4String ID;

    ID = "TED";
    (*store)[ID] = G4AttDef(ID,"Total Energy Deposit",
			    "Physics","G4BestUnit","G4double");
    ID = "PDS";
    (*store)[ID] = G4AttDef(ID,"Process Defined Step",
			    "Bookkeeping","","G4String");
    ID = "PTDS";
    (*store)[ID] = G4AttDef(ID,"Process Type Defined Step",
			    "Bookkeeping","","G4String");
  }
  return store;
}

std::vector<G4AttValue>* G4RichTrajectoryPoint::CreateAttValues() const
{
  std::vector<G4AttValue>* values = G4TrajectoryPoint::CreateAttValues();

  values->push_back(G4AttValue("TED",G4BestUnit(fTotEDep,"Energy"),""));

  if (fpProcess) {
    values->push_back
      (G4AttValue("PDS",fpProcess->GetProcessName(),""));
    values->push_back
      (G4AttValue
       ("PTDS",G4VProcess::GetProcessTypeName(fpProcess->GetProcessType()),
	""));
  } else {
    values->push_back
      (G4AttValue("PDS","User Defined Limit in Current Volume",""));
    values->push_back(G4AttValue("PTDS","User",""));
  }

#ifdef G4ATTDEBUG
  G4cout << G4AttCheck(values,GetAttDefs());
#endif

  return values;
}
