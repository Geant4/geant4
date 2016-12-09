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
// $Id: G4VScoringMesh.cc 99262 2016-09-09 13:18:19Z gcosmo $
//
// ---------------------------------------------------------------------
// Modifications                                                        
// 17-Apr-2012 T.Aso SetSize() and SetNumberOfSegments() is not allowed 
//                   to call twice in same geometrical mesh. Add warning
//                   message to notify.                                 
//                                                                      
// ---------------------------------------------------------------------

#include "G4VScoringMesh.hh"
#include "G4THitsMap.hh"
#include "G4SystemOfUnits.hh"
#include "G4VPhysicalVolume.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4VSDFilter.hh"
#include "G4SDManager.hh"

G4VScoringMesh::G4VScoringMesh(const G4String& wName)
  : fWorldName(wName),fCurrentPS(nullptr),fConstructed(false),fActive(true),fShape(undefinedMesh),
    fRotationMatrix(nullptr), fMFD(new G4MultiFunctionalDetector(wName)),
    verboseLevel(0),sizeIsSet(false),nMeshIsSet(false),
    fDrawUnit(""), fDrawUnitValue(1.), fMeshElementLogical(nullptr),
    fParallelWorldProcess(nullptr), fGeometryHasBeenDestroyed(false)
{
  G4SDManager::GetSDMpointer()->AddNewDetector(fMFD);

  fSize[0] = fSize[1] = fSize[2] = 0.;
  fNSegment[0] = fNSegment[1] = fNSegment[2] = 1;
  fDivisionAxisNames[0] = fDivisionAxisNames[1] = fDivisionAxisNames[2] = "";
}

G4VScoringMesh::~G4VScoringMesh() {
  ;
}

void G4VScoringMesh::ResetScore() {
  if(verboseLevel > 9) G4cout << "G4VScoringMesh::ResetScore() is called." << G4endl;
  for(auto mp : fMap)
  {
    if(verboseLevel > 9) G4cout << "G4VScoringMesh::ResetScore()" << mp.first << G4endl;
    mp.second->clear();
  }
}

void G4VScoringMesh::SetSize(G4double size[3]) {
  if ( !sizeIsSet ){
    for(int i = 0; i < 3; i++) fSize[i] = size[i];
    sizeIsSet = true;
  }else{
    G4String message = "   The size of scoring mesh can not be changed.";
    G4Exception("G4VScoringMesh::SetSize()",
                "DigiHitsUtilsScoreVScoringMesh000", JustWarning,
                message);
  }
}
G4ThreeVector G4VScoringMesh::GetSize() const {
  if(sizeIsSet)
    return G4ThreeVector(fSize[0], fSize[1], fSize[2]);
  else
    return G4ThreeVector(0., 0., 0.);
}
void G4VScoringMesh::SetCenterPosition(G4double centerPosition[3]) {
  fCenterPosition = G4ThreeVector(centerPosition[0], centerPosition[1], centerPosition[2]);
}
void G4VScoringMesh::SetNumberOfSegments(G4int nSegment[3]) {
  if ( !nMeshIsSet ){
    for(int i = 0; i < 3; i++) fNSegment[i] = nSegment[i];
    nMeshIsSet = true;
  } else {
    G4String message = "   The size of scoring segments can not be changed.";
    G4Exception("G4VScoringMesh::SetNumberOfSegments()",
                "DigiHitsUtilsScoreVScoringMesh000", JustWarning,
                message);
  }
}
void G4VScoringMesh::GetNumberOfSegments(G4int nSegment[3]) {
  for(int i = 0; i < 3; i++) nSegment[i] = fNSegment[i];
}
void G4VScoringMesh::RotateX(G4double delta) {
  if(!fRotationMatrix) fRotationMatrix = new G4RotationMatrix();
  fRotationMatrix->rotateX(delta);
}

void G4VScoringMesh::RotateY(G4double delta) {
  if(!fRotationMatrix) fRotationMatrix = new G4RotationMatrix();
  fRotationMatrix->rotateY(delta);
}

void G4VScoringMesh::RotateZ(G4double delta) {
  if(!fRotationMatrix) fRotationMatrix = new G4RotationMatrix();
  fRotationMatrix->rotateZ(delta);
}

void G4VScoringMesh::SetPrimitiveScorer(G4VPrimitiveScorer * prs) {

  if(!ReadyForQuantity())
  {
    G4cerr << "ERROR : G4VScoringMesh::SetPrimitiveScorer() : "
           << prs->GetName() 
           << " does not yet have mesh size or number of bins. Set them first." << G4endl
           << "This Method is ignored." << G4endl;
    return;
  }
  if(verboseLevel > 0) G4cout << "G4VScoringMesh::SetPrimitiveScorer() : "
			      << prs->GetName() << " is registered."
			      << " 3D size: ("
			      << fNSegment[0] << ", "
			      << fNSegment[1] << ", "
			      << fNSegment[2] << ")" << G4endl;

  prs->SetNijk(fNSegment[0], fNSegment[1], fNSegment[2]);
  fCurrentPS = prs;
  fMFD->RegisterPrimitive(prs);
  G4THitsMap<G4StatDouble> * map = new G4THitsMap<G4StatDouble>(fWorldName, prs->GetName());
  fMap[prs->GetName()] = map;
}

void G4VScoringMesh::SetFilter(G4VSDFilter * filter) {

  if(!fCurrentPS) {
    G4cerr << "ERROR : G4VScoringMesh::SetSDFilter() : a quantity must be defined first. This method is ignored." << G4endl;
    return;
  }
  if(verboseLevel > 0) G4cout << "G4VScoringMesh::SetFilter() : "
			      << filter->GetName()
			      << " is set to "
			      << fCurrentPS->GetName() << G4endl;

  G4VSDFilter* oldFilter = fCurrentPS->GetFilter();
  if(oldFilter)
  {
    G4cout << "WARNING : G4VScoringMesh::SetFilter() : " << oldFilter->GetName() 
           << " is overwritten by " << filter->GetName() << G4endl;
  }
  fCurrentPS->SetFilter(filter);
}

void G4VScoringMesh::SetCurrentPrimitiveScorer(const G4String & name) {
  fCurrentPS = GetPrimitiveScorer(name);
  if(!fCurrentPS) {
    G4cerr << "ERROR : G4VScoringMesh::SetCurrentPrimitiveScorer() : The primitive scorer <"
	   << name << "> does not found." << G4endl;
  }
}

G4bool G4VScoringMesh::FindPrimitiveScorer(const G4String & psname) {
  MeshScoreMap::iterator itr = fMap.find(psname);
  if(itr == fMap.end()) return false;
  return true;
}

G4String G4VScoringMesh::GetPSUnit(const G4String & psname) {
  MeshScoreMap::iterator itr = fMap.find(psname);
  if(itr == fMap.end()) {
    return G4String("");
  } else {
    return GetPrimitiveScorer(psname)->GetUnit();
  }
}

G4String G4VScoringMesh::GetCurrentPSUnit(){
  G4String unit = "";
  if(!fCurrentPS) {
      G4String msg = "ERROR : G4VScoringMesh::GetCurrentPSUnit() : ";
      msg += " Current primitive scorer is null.";
      G4cerr << msg << G4endl;
  }else{
     unit =  fCurrentPS->GetUnit();
  }
  return unit;
}

void  G4VScoringMesh::SetCurrentPSUnit(const G4String& unit){
  if(!fCurrentPS) {
      G4String msg = "ERROR : G4VScoringMesh::GetCurrentPSUnit() : ";
      msg += " Current primitive scorer is null.";
      G4cerr << msg << G4endl;
  }else{
      fCurrentPS->SetUnit(unit);
  }
}

G4double G4VScoringMesh::GetPSUnitValue(const G4String & psname) {
  MeshScoreMap::iterator itr = fMap.find(psname);
  if(itr == fMap.end()) {
    return 1.;
  } else {
    return GetPrimitiveScorer(psname)->GetUnitValue();
  }
}

void G4VScoringMesh::GetDivisionAxisNames(G4String divisionAxisNames[3]) {
  for(int i = 0; i < 3; i++) divisionAxisNames[i] = fDivisionAxisNames[i];
}

G4VPrimitiveScorer * G4VScoringMesh::GetPrimitiveScorer(const G4String & name) {
  if(!fMFD) return nullptr;

  G4int nps = fMFD->GetNumberOfPrimitives();
  for(G4int i = 0; i < nps; i++) {
    G4VPrimitiveScorer * prs = fMFD->GetPrimitive(i);
    if(name == prs->GetName()) return prs;
  }

  return nullptr;
}
void G4VScoringMesh::List() const {

  G4cout << " # of segments: ("
	 << fNSegment[0] << ", "
	 << fNSegment[1] << ", "
	 << fNSegment[2] << ")"
	 << G4endl;
  G4cout << " displacement: ("
	 << fCenterPosition.x()/cm << ", "
	 << fCenterPosition.y()/cm << ", "
	 << fCenterPosition.z()/cm << ") [cm]"
	 << G4endl;
  if(fRotationMatrix != 0) {
    G4cout << " rotation matrix: "
	   << fRotationMatrix->xx() << "  "
	   << fRotationMatrix->xy() << "  "
	   << fRotationMatrix->xz() << G4endl
	   << "                  "
	   << fRotationMatrix->yx() << "  "
	   << fRotationMatrix->yy() << "  "
	   << fRotationMatrix->yz() << G4endl
	   << "                  "
	   << fRotationMatrix->zx() << "  "
	   << fRotationMatrix->zy() << "  "
	   << fRotationMatrix->zz() << G4endl;
  }


  G4cout << " registered primitve scorers : " << G4endl;
  G4int nps = fMFD->GetNumberOfPrimitives();
  G4VPrimitiveScorer * prs;
  for(int i = 0; i < nps; i++) {
    prs = fMFD->GetPrimitive(i);
    G4cout << "   " << i << "  " << prs->GetName();
    if(prs->GetFilter() != 0) G4cout << "     with  " << prs->GetFilter()->GetName();
    G4cout << G4endl;
  }
}

void G4VScoringMesh::Dump() {
  G4cout << "scoring mesh name: " << fWorldName << G4endl;
  G4cout << "# of G4THitsMap : " << fMap.size() << G4endl;
  for(auto mp : fMap)
  {
    G4cout << "[" << mp.first << "]" << G4endl;
    mp.second->PrintAllHits();
  }
  G4cout << G4endl;

}


void G4VScoringMesh::DrawMesh(const G4String& psName,G4VScoreColorMap* colorMap,G4int axflg)
{
  fDrawPSName = psName;
  MeshScoreMap::const_iterator fMapItr = fMap.find(psName);
  if(fMapItr!=fMap.end()) {
    fDrawUnit = GetPSUnit(psName);
    fDrawUnitValue = GetPSUnitValue(psName);
    Draw(fMapItr->second, colorMap,axflg);
  } else {
    G4cerr << "Scorer <" << psName << "> is not defined. Method ignored." << G4endl;
  }
}

void G4VScoringMesh::DrawMesh(const G4String& psName,G4int idxPlane,G4int iColumn,G4VScoreColorMap* colorMap)
{
  fDrawPSName = psName;
  MeshScoreMap::const_iterator fMapItr = fMap.find(psName);
  if(fMapItr!=fMap.end()) {
    fDrawUnit = GetPSUnit(psName);
    fDrawUnitValue = GetPSUnitValue(psName);
    DrawColumn(fMapItr->second,colorMap,idxPlane,iColumn);
  } else {
    G4cerr << "Scorer <" << psName << "> is not defined. Method ignored." << G4endl;
  }
}

void G4VScoringMesh::Accumulate(G4THitsMap<G4double> * map)
{
  G4String psName = map->GetName();
  MeshScoreMap::const_iterator fMapItr = fMap.find(psName);
  *(fMapItr->second) += *map;

  if(verboseLevel > 9) {
    G4cout << G4endl;
    G4cout << "G4VScoringMesh::Accumulate()" << G4endl;
    G4cout << "  PS name : " << psName << G4endl;
    if(fMapItr == fMap.end()) {
      G4cout << "  " << psName << " was not found." << G4endl;
    } else {
      G4cout << "  map size : " << map->GetSize() << G4endl;
      map->PrintAllHits();
    }
    G4cout << G4endl;
  }
}

void G4VScoringMesh::Accumulate(G4THitsMap<G4StatDouble> * map)
{
  G4String psName = map->GetName();
  MeshScoreMap::const_iterator fMapItr = fMap.find(psName);
  *(fMapItr->second) += *map;

  if(verboseLevel > 9) {
    G4cout << G4endl;
    G4cout << "G4VScoringMesh::Accumulate()" << G4endl;
    G4cout << "  PS name : " << psName << G4endl;
    if(fMapItr == fMap.end()) {
      G4cout << "  " << psName << " was not found." << G4endl;
    } else {
      G4cout << "  map size : " << map->GetSize() << G4endl;
      map->PrintAllHits();
    }
    G4cout << G4endl;
  }
}

void G4VScoringMesh::Construct(G4VPhysicalVolume* fWorldPhys)
{
  if(fConstructed) {
    if(fGeometryHasBeenDestroyed) {
      SetupGeometry(fWorldPhys);
      fGeometryHasBeenDestroyed = false;
    }
    if(verboseLevel > 0)
      G4cout << fWorldPhys->GetName() << " --- All quantities are reset."
      << G4endl;
    ResetScore();
  }
  else {
    fConstructed = true;
    SetupGeometry(fWorldPhys);
  }
}

void G4VScoringMesh::WorkerConstruct(G4VPhysicalVolume* fWorldPhys)
{
  if(fConstructed) {
    if(fGeometryHasBeenDestroyed) {
      fMeshElementLogical->SetSensitiveDetector(fMFD);
      fGeometryHasBeenDestroyed = false;
    }

    if(verboseLevel > 0)
      G4cout << fWorldPhys->GetName() << " --- All quantities are reset." << G4endl;
    ResetScore();

  } else {
    fConstructed = true;
    fMeshElementLogical->SetSensitiveDetector(fMFD);
  }
}

void G4VScoringMesh::Merge(const G4VScoringMesh * scMesh)
{
  const MeshScoreMap scMap = scMesh->GetScoreMap();

  MeshScoreMap::const_iterator fMapItr = fMap.begin();
  MeshScoreMap::const_iterator mapItr = scMap.begin();
  for(; fMapItr != fMap.end(); fMapItr++) {
    if(verboseLevel > 9) G4cout << "G4VScoringMesh::Merge()" << fMapItr->first << G4endl;
    *(fMapItr->second) += *(mapItr->second);
    mapItr++;
  }
}

