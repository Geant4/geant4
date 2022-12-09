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
#ifndef G4ScoringRealWorld_h
#define G4ScoringRealWorld_h 1

#include "globals.hh"
#include "G4VScoringMesh.hh"
class G4VPhysicalVolume;
class G4LogicalVolume;

class G4ScoringRealWorld : public G4VScoringMesh
{
 public:
  G4ScoringRealWorld(G4String lvName);
  ~G4ScoringRealWorld() override = default;

 public:
  void List() const override;

  //++++++++++ visualization method not yet implemented
  void Draw(RunScore* /*map*/, G4VScoreColorMap* /*colorMap*/,
                    G4int /*axflg=111*/) override
  {}

  void DrawColumn(RunScore* /*map*/, G4VScoreColorMap* /*colorMap*/,
                          G4int /*idxProj*/, G4int /*idxColumn*/) override
  {}

 protected:
  // construct this mesh
  void SetupGeometry(G4VPhysicalVolume*) override;

 protected:
  G4String logVolName;
};

#endif
