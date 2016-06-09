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
// $Id: G4RichTrajectoryPoint.cc,v 1.3 2006/09/27 20:42:52 asaim Exp $
// GEANT4 tag $Name: geant4-09-01 $
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
  fpAuxiliaryPointVector(0),
  fTotEDep(0.),
  fpProcess(0),
  fPreStepPointGlobalTime(0),
  fPostStepPointGlobalTime(0)
{}

G4RichTrajectoryPoint::G4RichTrajectoryPoint(const G4Track* aTrack):
  G4TrajectoryPoint(aTrack->GetPosition()),
  fpAuxiliaryPointVector(0),
  fTotEDep(0.),
  fpProcess(0),
  fPreStepPointGlobalTime(aTrack->GetGlobalTime()),
  fPostStepPointGlobalTime(aTrack->GetGlobalTime())
{}

G4RichTrajectoryPoint::G4RichTrajectoryPoint(const G4Step* aStep):
  G4TrajectoryPoint(aStep->GetPostStepPoint()->GetPosition()),
  fpAuxiliaryPointVector(aStep->GetPointerToVectorOfAuxiliaryPoints()),
  fTotEDep(aStep->GetTotalEnergyDeposit())
{
  G4StepPoint* preStepPoint = aStep->GetPreStepPoint();
  G4StepPoint* postStepPoint = aStep->GetPostStepPoint();
  fpProcess = postStepPoint->GetProcessDefinedStep();
  fPreStepPointGlobalTime = preStepPoint->GetGlobalTime();
  fPostStepPointGlobalTime = postStepPoint->GetGlobalTime();

  /*
  G4cout << "fpAuxiliaryPointVector "
	 << (void*) fpAuxiliaryPointVector;
  G4cout << ": ";
  if (fpAuxiliaryPointVector) {
    G4cout << "size: " << fpAuxiliaryPointVector->size();
    for (size_t i = 0; i < fpAuxiliaryPointVector->size(); ++i)
      G4cout << "\n  " << (*fpAuxiliaryPointVector)[i];
  } else {
    G4cout << "non-existent";
  }
  G4cout << G4endl;

  static const G4Step* lastStep = 0;
  if (aStep && aStep == lastStep) {
    G4cout << "********* aStep is same as last" << G4endl;
  }
  lastStep = aStep;

  static std::vector<G4ThreeVector>*  lastAuxiliaryPointVector = 0;
  if (fpAuxiliaryPointVector &&
      fpAuxiliaryPointVector == lastAuxiliaryPointVector) {
    G4cout << "********* fpAuxiliaryPointVector is same as last" << G4endl;
  }
  lastAuxiliaryPointVector = fpAuxiliaryPointVector;
  */
}

G4RichTrajectoryPoint::G4RichTrajectoryPoint
(const G4RichTrajectoryPoint &right):
  G4TrajectoryPoint(right),
  fpAuxiliaryPointVector(right.fpAuxiliaryPointVector),
  fTotEDep(right.fTotEDep),
  fpProcess(right.fpProcess),
  fPreStepPointGlobalTime(right.fPreStepPointGlobalTime),
  fPostStepPointGlobalTime(right.fPostStepPointGlobalTime)
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
    ID = "PDS";
    (*store)[ID] = G4AttDef(ID,"Process Defined Step",
			    "Physics","","G4String");
    ID = "PTDS";
    (*store)[ID] = G4AttDef(ID,"Process Type Defined Step",
			    "Physics","","G4String");
    ID = "PreT";
    (*store)[ID] = G4AttDef(ID,"Pre-step-point global time",
			    "Physics","G4BestUnit","G4double");
    ID = "PostT";
    (*store)[ID] = G4AttDef(ID,"Post-step-point global time",
			    "Physics","G4BestUnit","G4double");
  }
  return store;
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

  values->push_back
    (G4AttValue("PreT",G4BestUnit(fPreStepPointGlobalTime,"Time"),""));

  values->push_back
    (G4AttValue("PostT",G4BestUnit(fPostStepPointGlobalTime,"Time"),""));

#ifdef G4ATTDEBUG
  G4cout << G4AttCheck(values,GetAttDefs());
#endif

  return values;
}
