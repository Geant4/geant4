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
// $Id: G4VScoringMesh.cc,v 1.5 2007-08-28 05:26:56 akimura Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "G4VScoringMesh.hh"
#include "G4VPhysicalVolume.hh"
#include "G4MultiFunctionalDetector.hh"

G4VScoringMesh::G4VScoringMesh(G4String wName)
  : fWorldName(wName),fConstructed(false),fActive(true),
    fRotationMatrix(NULL), fMFD(new G4MultiFunctionalDetector(wName)),
    fScoringMeshName(wName) 
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

#include "G4PSEnergyDeposit.hh"
#include "G4PSDoseDeposit.hh"

void G4VScoringMesh::CreatePrimitiveScorer(PS psType, G4String & name) {

  G4VPrimitiveScorer * ps;
  switch(psType) {
  case EnergyDeposit:
    ps = new G4PSEnergyDeposit(name);
    break;
  case DoseDeposit: 
    ps = new G4PSDoseDeposit(name);
    break;
  default:
    return;
  }
  fMFD->RegisterPrimitive(ps);
  
}

#include "G4SDParticleFilter.hh"
#include "G4SDChargedFilter.hh"

void G4VScoringMesh::CreateSDFilter(G4String & psName,
				    FILTER filterType,
				    G4String & filterName,
				    std::vector<G4String> & parameter) {

  G4VSDFilter * filter;
  switch(filterType) {
  case Particle:
    filter = new G4SDParticleFilter(filterName, parameter);
  case Charged:
    filter = new G4SDChargedFilter(filterName);
  default:
    return;
  }

  G4VPrimitiveScorer * ps = GetPrimitiveScorer(psName);
  if(ps == NULL) return;

  ps->SetFilter(filter);
  
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

