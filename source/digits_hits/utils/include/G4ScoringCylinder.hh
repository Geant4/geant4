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
//

#ifndef G4ScoringCylinder_h
#define G4ScoringCylinder_h 1

#include "globals.hh"
#include "G4VScoringMesh.hh"
class G4VPhysicalVolume;
class G4LogicalVolume;
class G4VPrimitiveScorer;

#include <vector>

class G4ScoringCylinder : public G4VScoringMesh
{
 public:
  G4ScoringCylinder(G4String wName);
  ~G4ScoringCylinder() override = default;


 public:
  void List() const override;
  void Draw(RunScore* map, G4VScoreColorMap* colorMap,
                    G4int axflg = 111) override;
  void DrawColumn(RunScore* map, G4VScoreColorMap* colorMap,
                          G4int idxProj, G4int idxColumn) override;

  void SetRMin(G4double rMin) { fSize[0] = rMin; }
  void SetRMax(G4double rMax) { fSize[1] = rMax; }
  void SetZSize(G4double zSize) { fSize[2] = zSize; }  // half height

  void RegisterPrimitives(std::vector<G4VPrimitiveScorer*>& vps);

  // get 3D index (z,phi,r) from sequential index
  void GetRZPhi(G4int index, G4int q[3]) const;

 protected:
  void SetupGeometry(G4VPhysicalVolume* fWorldPhys) override;

 private:
  // Xin Dong 09302011 for Scorers
  enum IDX
  {
    IZ,
    IPHI,
    IR
  };

  void DumpVolumes();
  void DumpSolids(G4int);
  void DumpLogVols(G4int);
  void DumpPhysVols(G4int);
};

#endif
