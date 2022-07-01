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

#include "globals.hh"
#include "G4DNAMakeReaction.hh"
#include <G4DNAMolecularReactionTable.hh>
#include <memory>
#include "G4DNAIndependentReactionTimeModel.hh"
#include "G4DNAIndependentReactionTimeStepper.hh"
#include "G4DiffusionControlledReactionModel.hh"
#include "G4Molecule.hh"
#include "G4VDNAReactionModel.hh"
#include "G4ChemicalMoleculeFinder.hh"

G4DNAIndependentReactionTimeModel::G4DNAIndependentReactionTimeModel(
  const G4String& name)
  : G4DNAIndependentReactionTimeModel(
      name, std::make_unique<G4DNAIndependentReactionTimeStepper>(),
      std::make_unique<G4DNAMakeReaction>())
{}

G4DNAIndependentReactionTimeModel::G4DNAIndependentReactionTimeModel(
  const G4String& name, std::unique_ptr<G4VITTimeStepComputer> pTimeStepper,
  std::unique_ptr<G4VITReactionProcess> pReactionProcess)
  : G4VITStepModel(std::move(pTimeStepper), std::move(pReactionProcess), name)
{
  fType1 = G4Molecule::ITType();
  fType2 = G4Molecule::ITType();
}

G4DNAIndependentReactionTimeModel::~G4DNAIndependentReactionTimeModel() =
  default;

void G4DNAIndependentReactionTimeModel::Initialize()
{
  if(fpReactionTable == nullptr)
  {
    SetReactionTable(G4DNAMolecularReactionTable::GetReactionTable());
  }

  if(!fpReactionModel)
  {
    fpReactionModel = std::make_unique<G4DiffusionControlledReactionModel>();
  }

  fpReactionModel->SetReactionTable(
    (const G4DNAMolecularReactionTable*) fpReactionTable);

  ((G4DNAMakeReaction*) fpReactionProcess.get())
    ->SetReactionModel(fpReactionModel.get());

  ((G4DNAMakeReaction*) fpReactionProcess.get())
    ->SetTimeStepComputer(fpTimeStepper.get());

  ((G4DNAIndependentReactionTimeStepper*) fpTimeStepper.get())
    ->SetReactionModel(fpReactionModel.get());

  ((G4DNAIndependentReactionTimeStepper*) fpTimeStepper.get())
    ->SetReactionProcess((fpReactionProcess).get());

  G4ChemicalMoleculeFinder::Instance()->Clear();
  G4ChemicalMoleculeFinder::Instance()->SetOctreeUsed(true);
  G4VITStepModel::Initialize();
}
void G4DNAIndependentReactionTimeModel::PrintInfo()
{
#ifdef G4VERBOSE
  if(G4Threading::IsMultithreadedApplication())
  {
    if(G4Threading::G4GetThreadId() == 0)
    {
      G4VITStepModel::PrintInfo();
      G4cout << G4endl;
      G4cout << fName << " will be used ==========================" << G4endl;
      G4cout << G4endl;
      G4cout << "==============================="
                "========================================"
             << G4endl;
      G4cout << G4endl;
    }
  }
  else
  {
    G4VITStepModel::PrintInfo();
    G4cout << G4endl;
    G4cout << fName << " will be used ==========================" << G4endl;
    G4cout << G4endl;
    G4cout << "==============================="
              "========================================"
           << G4endl;
    G4cout << G4endl;
  }
#endif
}

void G4DNAIndependentReactionTimeModel::SetReactionModel(
  G4VDNAReactionModel* pReactionModel)
{
  fpReactionModel.reset(pReactionModel);
}

G4VDNAReactionModel* G4DNAIndependentReactionTimeModel::GetReactionModel()
{
  return fpReactionModel.get();
}