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
// $Id$
//

#ifndef G4VScoringMesh_h
#define G4VScoringMesh_h 1

#include "globals.hh"
#include "G4THitsMap.hh"
#include "G4RotationMatrix.hh"

class G4VPhysicalVolume;
class G4MultiFunctionalDetector;
class G4VPrimitiveScorer;
class G4VSDFilter;
class G4VScoreColorMap;

#include <map>

enum MeshShape { boxMesh, cylinderMesh, sphereMesh };
typedef std::map<G4String,G4THitsMap<G4double>* > MeshScoreMap;
// class description:
//
//  This class represents a parallel world for interactive scoring purposes.
//

class G4VScoringMesh
{
  public:
  G4VScoringMesh(const G4String& wName);
  virtual ~G4VScoringMesh();

  public: // with description
  // a pure virtual function to construct this mesh geometry
  virtual void Construct(G4VPhysicalVolume* fWorldPhys)=0;
  // list infomration of this mesh 
  virtual void List() const;
  
  public: // with description
  // get the world name
  inline const G4String& GetWorldName() const
  { return fWorldName; }
  // get whether this mesh is active or not
  inline G4bool IsActive() const
  { return fActive; }
  // set an activity of this mesh
  inline void Activate(G4bool vl = true)
  { fActive = vl; }
  // get the shape of this mesh
  inline MeshShape GetShape() const
  { return fShape; }
  // accumulate hits in a registered primitive scorer
  inline void Accumulate(G4THitsMap<G4double> * map);
  // dump information of primitive socrers registered in this mesh
  void Dump();
  // draw a projected quantity on a current viewer
  void DrawMesh(const G4String& psName,G4VScoreColorMap* colorMap,G4int axflg=111);
  // draw a column of a quantity on a current viewer
  void DrawMesh(const G4String& psName,G4int idxPlane,G4int iColumn,G4VScoreColorMap* colorMap);
  // draw a projected quantity on a current viewer
  virtual void Draw(std::map<G4int, G4double*> * map, G4VScoreColorMap* colorMap, G4int axflg=111) = 0;
  // draw a column of a quantity on a current viewer
  virtual void DrawColumn(std::map<G4int, G4double*> * map, G4VScoreColorMap* colorMap,
			  G4int idxProj, G4int idxColumn) = 0;
  // reset registered primitive scorers
  void ResetScore();

  // set size of this mesh
  void SetSize(G4double size[3]);
  // get size of this mesh
  G4ThreeVector GetSize() const;
  // set position of center of this mesh
  void SetCenterPosition(G4double centerPosition[3]);
  // get position of center of this mesh
  G4ThreeVector GetTranslation() const {return fCenterPosition;}
  // set a rotation angle around the x axis
  void RotateX(G4double delta);
  // set a rotation angle around the y axis
  void RotateY(G4double delta);
  // set a rotation angle around the z axis
  void RotateZ(G4double delta);
  // get a rotation matrix
  G4RotationMatrix GetRotationMatrix() const {
    if(fRotationMatrix) return *fRotationMatrix;
    else return G4RotationMatrix::IDENTITY;
  }
  // set number of segments of this mesh
  void SetNumberOfSegments(G4int nSegment[3]);
  // get number of segments of this mesh
  void GetNumberOfSegments(G4int nSegment[3]);

  // register a primitive scorer to the MFD & set it to the current primitive scorer
  void SetPrimitiveScorer(G4VPrimitiveScorer * ps);
  // register a filter to a current primtive scorer
  void SetFilter(G4VSDFilter * filter);
  // set a primitive scorer to the current one by the name
  void SetCurrentPrimitiveScorer(const G4String & name);
  // find registered primitive scorer by the name
  G4bool FindPrimitiveScorer(const G4String & psname);
  // get whether current primitive scorer is set or not
  G4bool IsCurrentPrimitiveScorerNull() {
    if(fCurrentPS == NULL) return true;
    else return false;
  }
  // get unit of primitive scorer by the name
  G4String GetPSUnit(const G4String & psname);
  // get unit of current primitive scorer
  G4String GetCurrentPSUnit();
  // set unit of current primitive scorer
  void SetCurrentPSUnit(const G4String& unit);
  // get unit value of primitive scorer by the name
  G4double GetPSUnitValue(const G4String & psname);
  // set PS name to be drawn
  void SetDrawPSName(const G4String & psname) {fDrawPSName = psname;}

  // get axis names of the hierarchical division in the divided order
  void GetDivisionAxisNames(G4String divisionAxisNames[3]);

  // set current  primitive scorer to NULL
  void SetNullToCurrentPrimitiveScorer() {fCurrentPS = NULL;}
  // set verbose level
  inline void SetVerboseLevel(G4int vl) 
  { verboseLevel = vl; }
  // get the primitive scorer map
  MeshScoreMap GetScoreMap() {return fMap;}
  // get whether this mesh setup has been ready
  inline G4bool ReadyForQuantity() const
  { return (sizeIsSet && nMeshIsSet); }

protected:
  // get registered primitive socrer by the name
  G4VPrimitiveScorer * GetPrimitiveScorer(const G4String & name);

protected:
  G4String  fWorldName;
  G4VPrimitiveScorer * fCurrentPS;
  G4bool    fConstructed;
  G4bool    fActive;
  MeshShape fShape;

  G4double fSize[3];
  G4ThreeVector fCenterPosition;
  G4RotationMatrix * fRotationMatrix;
  G4int fNSegment[3];

  std::map<G4String, G4THitsMap<G4double>* > fMap;
  G4MultiFunctionalDetector * fMFD;

  G4int verboseLevel;

  G4bool sizeIsSet;
  G4bool nMeshIsSet;

  G4String fDrawUnit;
  G4double fDrawUnitValue;
  G4String fDrawPSName;

  G4String fDivisionAxisNames[3];
};

void G4VScoringMesh::Accumulate(G4THitsMap<G4double> * map)
{
  G4String psName = map->GetName();
  std::map<G4String, G4THitsMap<G4double>* >::const_iterator fMapItr = fMap.find(psName);
  *(fMapItr->second) += *map;

  if(verboseLevel > 9) {
    G4cout << G4endl;
    G4cout << "G4VScoringMesh::Accumulate()" << G4endl;
    G4cout << "  PS name : " << psName << G4endl;
    if(fMapItr == fMap.end()) {
      G4cout << "  "
	     << psName << " was not found." << G4endl;
    } else {
      G4cout << "  map size : " << map->GetSize() << G4endl;
      map->PrintAllHits();
    }
    G4cout << G4endl;
  }
}

#endif

