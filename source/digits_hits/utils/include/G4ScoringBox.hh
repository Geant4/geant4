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
// $Id: G4ScoringBox.hh,v 1.3 2007-08-10 08:36:47 akimura Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef G4ScoringBox_h
#define G4ScoringBox_h 1

#include "globals.hh"
#include "G4VScoringMesh.hh"
#include "G4RotationMatrix.hh"
class G4VPhysicalVolume;
class G4LogicalVolume;
class G4VPrimitiveScorer;

#include <vector>

class G4ScoringBox : public G4VScoringMesh
{
  public:
      G4ScoringBox(G4String wName);
      ~G4ScoringBox();

  public:
      virtual void Construct(G4VPhysicalVolume* fWorldPhys);
      virtual void List() const;

  void SetSize(G4double size[3]) {
    for(int i = 0; i < 3; i++) fSize[i] = size[i];
  }
  void SetXSize(G4double xSize) {fSize[0] = xSize;}
  void SetYSize(G4double ySize) {fSize[1] = ySize;}
  void SetZSize(G4double zSize) {fSize[2] = zSize;}

  void SetCenterPosition(G4double centerPosition[3]) {
    for(int i = 0; i < 3; i++) fCenterPosition[i] = centerPosition[i];
  }
  void SetXCenterPosition(G4double xCenterPosition) {fCenterPosition[0] = xCenterPosition;}
  void SetYCenterPosition(G4double yCenterPosition) {fCenterPosition[1] = yCenterPosition;}
  void SetZCenterPosition(G4double zCenterPosition) {fCenterPosition[2] = zCenterPosition;}

  void SetNumberOfSegment(G4int nSegment[3]) {
    for(int i = 0; i < 3; i++) fNSegment[i] = nSegment[i];
}
  void SetNumberOfXSegment(G4int nXSegment) {fNSegment[0] = nXSegment;}
  void SetNumberOfYSegment(G4int nYSegment) {fNSegment[1] = nYSegment;}
  void SetNumberOfZSegment(G4int nZSegment) {fNSegment[2] = nZSegment;}

  void SetSegmentDirection(G4int dir) {fSegmentDirection = dir;}
  void SetRotationMatrix(G4RotationMatrix * rmat) {fRotationMatrix = rmat;}
  void SetSegmentPositions(std::vector<G4double> & sp) {fSegmentPositions = sp;}
  void RegisterPrimitives(std::vector<G4VPrimitiveScorer *> & vps);

private:
  G4double fSize[3];
  G4double fCenterPosition[3];
  G4RotationMatrix * fRotationMatrix;
  G4int fNSegment[3];
  G4int fSegmentDirection; // 0: x, 1: y, 2: z
  std::vector<G4double> fSegmentPositions;
  G4LogicalVolume * fMeshElementLogical;
  
  void SetupGeometry(G4VPhysicalVolume * fWorldPhys);
  void GetSegmentOrder(G4int segDir, G4int nseg[3], G4int segOrd[3], G4double segfact[3][3]);
};




#endif

