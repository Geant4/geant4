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
// $Id: G4VScoringMesh.cc,v 1.13 2007-08-29 01:42:17 akimura Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "G4VScoringMesh.hh"
#include "G4VPhysicalVolume.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4VSDFilter.hh"

G4VScoringMesh::G4VScoringMesh(G4String wName)
  : fWorldName(wName),fConstructed(false),fActive(true),
    fRotationMatrix(NULL), fMFD(new G4MultiFunctionalDetector(wName)),
    verboseLevel(0)
{
  fSize[0] = fSize[1] = fSize[2] = 0.*cm;
  fCenterPosition[0] = fCenterPosition[1] = fCenterPosition[2] = 0.*cm;
  fNSegment[0] = fNSegment[1] = fNSegment[2] = 1;
}

G4VScoringMesh::~G4VScoringMesh()
{
  //delete fMFD;
}

//void G4VScoringMesh::Construct(G4VPhysicalVolume* fWorldPhys)
//{
//  if(fConstructed) 
//  {
//    G4cerr << fWorldPhys->GetName() << G4endl;
//    G4Exception(fWorldName+" has already been built.");
//  }
//  fConstructed = true;
//
//}

void G4VScoringMesh::SetSize(G4double size[3]) {
  for(int i = 0; i < 3; i++) fSize[i] = size[i];
}
void G4VScoringMesh::SetCenterPosition(G4double centerPosition[3]) {
  for(int i = 0; i < 3; i++) fCenterPosition[i] = centerPosition[i];
}
void G4VScoringMesh::SetNumberOfSegments(G4int nSegment[3]) {
  for(int i = 0; i < 3; i++) fNSegment[i] = nSegment[i];
}

void G4VScoringMesh::RotateX(G4double delta) {
  if(fRotationMatrix == NULL) fRotationMatrix = new G4RotationMatrix();
  fRotationMatrix->rotateX(delta);
}

void G4VScoringMesh::RotateY(G4double delta) {
  if(fRotationMatrix == NULL) fRotationMatrix = new G4RotationMatrix();
  fRotationMatrix->rotateY(delta);
}

void G4VScoringMesh::RotateZ(G4double delta) {
  if(fRotationMatrix == NULL) fRotationMatrix = new G4RotationMatrix();
  fRotationMatrix->rotateZ(delta);
}

void G4VScoringMesh::SetPrimitiveScorer(G4VPrimitiveScorer * ps) {

  ps->SetNijk(fNSegment[0], fNSegment[1], fNSegment[2]);
  fCurrentPS = ps;
  fMFD->RegisterPrimitive(ps);
  G4THitsMap<G4double> * map = new G4THitsMap<G4double>(fWorldName, ps->GetName());
  fMap[ps->GetName()] = map;
  if(verboseLevel > 9) G4cout << fMFD->GetNumberOfPrimitives() << G4endl;;
}

void G4VScoringMesh::SetFilter(G4VSDFilter * filter) {

  if(fCurrentPS == NULL) {
    //G4Exception("G4VScoringMesh::SetSDFilter() : Current primitive scorer has not been set.");
  }

  fCurrentPS->SetFilter(filter);
}

void G4VScoringMesh::SetCurrentPrimitiveScorer(G4String & name) {
  fCurrentPS = GetPrimitiveScorer(name);
  if(fCurrentPS == NULL) {
    //G4Execption("G4VScoringMesh::SetCurrentPrimitiveScorer() :  The primitive scorer ", name, " was not exist.");
  }
}

G4VPrimitiveScorer * G4VScoringMesh::GetPrimitiveScorer(G4String & name) {
  if(fMFD == NULL) return NULL;

  G4int nps = fMFD->GetNumberOfPrimitives();
  for(G4int i = 0; i < nps; i++) {
    G4VPrimitiveScorer * ps = fMFD->GetPrimitive(i);
    if(name == ps->GetName()) return ps;
  }

  return NULL;
}

void G4VScoringMesh::List() const {
  G4cout << " Size: " << G4endl;
  G4cout << " # of segments: " << G4endl;
  G4cout << " displacement: " << G4endl;
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

void G4VScoringMesh::Dump() {
  G4cout << "scoring mesh name: " << fWorldName << G4endl;

  G4cout << "# of G4THitsMap : " << fMap.size() << G4endl;
  std::map<G4String, G4THitsMap<G4double>* >::iterator itr = fMap.begin();
  for(; itr != fMap.end(); itr++) {
    G4cout << "[" << itr->first << "]" << G4endl;
    std::map<G4int, double*> * map = itr->second->GetMap();
    std::map<G4int, double*>::iterator itrMap  = map->begin();
    for(; itrMap != map->end(); itrMap++) {
      G4cout << itrMap->second << ", ";
    }
    G4cout << G4endl << G4endl;
  }
}

