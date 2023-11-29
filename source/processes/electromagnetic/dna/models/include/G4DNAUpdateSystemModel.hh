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
#ifndef G4DNAUPDATESYSTEMMODEL_HH
#define G4DNAUPDATESYSTEMMODEL_HH
#include "G4DNAMesh.hh"
#include "G4VUpdateSystemModel.hh"

class G4MolecularConfiguration;
class G4DNAMolecularReactionTable;
class G4DNAMolecularReactionData;

class G4DNAUpdateSystemModel : public G4VUpdateSystemModel
{
 public:
  using Index        = G4VDNAMesh::Index;
  using MolType      = const G4MolecularConfiguration*;
  using JumpingData  = std::pair<MolType, Index>;
  using ReactionData = const G4DNAMolecularReactionData;

  G4DNAUpdateSystemModel();
  ~G4DNAUpdateSystemModel() override = default;
  void UpdateSystem(const Index& index, const ReactionData& data);
  void UpdateSystem(const Index& index, const JumpingData& data);
  void SetMesh(G4DNAMesh*);
  void SetGlobalTime(const G4double& globalTime) { fGlobalTime = globalTime; }
  void SetVerbose(G4int verbose) { fVerbose = verbose; }

 private:
  void KillMolecule(const Index& index, MolType type);
  void CreateMolecule(const Index& index, MolType);
  void JumpTo(const Index& index, MolType type);
  void JumpIn(const Index& index, MolType);
  G4DNAMesh* fpMesh  = nullptr;
  G4int fVerbose = 0;
  G4double fGlobalTime = DBL_MAX;
};

#endif  // G4DNAUPDATESYSTEMMODEL_HH
