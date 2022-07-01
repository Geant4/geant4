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
// 20/2/2019
// Author : HoangTRAN

#ifndef G4DNAIndependentReactionTimeModel_hh
#define G4DNAIndependentReactionTimeModel_hh 1
#include "G4String.hh"
#include "G4VITStepModel.hh"
#include <G4ReferenceCast.hh>

class G4DNAMolecularReactionTable;
class G4VDNAReactionModel;
class G4DNAIndependentReactionTimeModel : public G4VITStepModel
{
 public:
  explicit G4DNAIndependentReactionTimeModel(
    const G4String& name = "DNAIndependentReactionTimeModel");
  G4DNAIndependentReactionTimeModel(
    const G4String& name, std::unique_ptr<G4VITTimeStepComputer> pTimeStepper,
    std::unique_ptr<G4VITReactionProcess> pReactionProcess);
  G4DNAIndependentReactionTimeModel(const G4DNAIndependentReactionTimeModel&) =
    delete;
  ~G4DNAIndependentReactionTimeModel() override;
  void PrintInfo() override;
  void Initialize() override;
  void SetReactionModel(G4VDNAReactionModel*);
  G4VDNAReactionModel* GetReactionModel();

 protected:
  const G4DNAMolecularReactionTable*& fMolecularReactionTable =
    reference_cast<const G4DNAMolecularReactionTable*>(fpReactionTable);
  std::unique_ptr<G4VDNAReactionModel> fpReactionModel;
};
#endif