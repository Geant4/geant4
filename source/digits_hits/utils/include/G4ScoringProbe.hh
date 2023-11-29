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

#ifndef G4ScoringProbe_h
#define G4ScoringProbe_h 1

#include "globals.hh"
#include "G4VScoringMesh.hh"

#include <vector>

class G4VPhysicalVolume;
class G4Material;

class G4ScoringProbe : public G4VScoringMesh
{
 public:
  G4ScoringProbe(G4String lvName, G4double half_size,
                 G4bool checkOverlap = false);
  ~G4ScoringProbe() override = default;

 public:
  void LocateProbe(G4ThreeVector pos)
  {
    posVec.push_back(pos);
    G4int nbin[] = { static_cast<G4int>(posVec.size()), 1, 1 };
    SetNumberOfSegments(nbin);
  }
  G4int GetNumberOfProbes() const { return (G4int)posVec.size(); }
  void SetProbeSize(G4double val) { probeSize = val; }
  G4double GetProbeSize() const { return probeSize; }
  G4bool SetMaterial(G4String val);

  void List() const override;

  //++++++++++ visualization method not yet implemented
  void Draw(RunScore* /*map*/, G4VScoreColorMap* /*colorMap*/,
                    G4int /*axflg=111*/) override
  {}
  void DrawColumn(RunScore* /*map*/, G4VScoreColorMap* /*colorMap*/,
                          G4int /*idxProj*/, G4int /*idxColumn*/) override
  {}

 protected:
  // construct scoring volume
  void SetupGeometry(G4VPhysicalVolume*) override;

 protected:
  G4String logVolName;
  std::vector<G4ThreeVector> posVec;
  G4double probeSize;
  G4bool chkOverlap;
  G4String layeredMaterialName;
  G4Material* layeredMaterial;
  G4String regName;
};

#endif
