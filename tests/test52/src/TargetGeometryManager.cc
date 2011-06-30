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

#include "TargetGeometryManager.hh"
#include "TargetComponent.hh"
#include "SemiInfiniteTarget.hh"
#include "AddOnTargetLayer.hh"


TargetGeometryManager::TargetGeometryManager() {

  layerRadius = 0.0;
  frontLayerThickness = 1.0 * cm;
  frontLayerMaterial = "";
  maxStepSize = 0.01 * mm;

  targetComponents = new Components;
}


TargetGeometryManager::~TargetGeometryManager() {

  delete targetComponents;
}


void TargetGeometryManager::Attach(TargetComponent* comp) {

  targetComponents -> push_back(comp);

  // Notify();
}


void TargetGeometryManager::Notify() {

  Components::reverse_iterator riter = targetComponents -> rbegin();
  Components::reverse_iterator riter_end = targetComponents -> rend();

  for(;riter != riter_end;riter++) {
     (*riter) -> GeometryUpdate(this);
  }
}


TargetComponent* TargetGeometryManager::GetFrontLayer() {

  Components::reverse_iterator iter = targetComponents -> rbegin();
  return (*iter);
}


G4double TargetGeometryManager::GetThickness(TargetComponent* comp) {
  
  Components::reverse_iterator iter = targetComponents -> rbegin();
  if((*iter) == comp) return frontLayerThickness;

  return 0.0;
}


G4String TargetGeometryManager::GetMaterial(TargetComponent* comp) {
  
  Components::reverse_iterator iter = targetComponents -> rbegin();
  if((*iter) == comp) return frontLayerMaterial;

  return "";
}


G4double TargetGeometryManager::GetPosition(TargetComponent* comp) {
  
  Components::reverse_iterator riter = targetComponents -> rbegin();
  Components::reverse_iterator riter_end = targetComponents -> rend();

  G4double pos = 0.0;
 
  for(;riter != riter_end;riter++) {
     if((*riter) == comp) {
        pos += 0.5 * (*riter) -> Thickness(); 
        break;
     }
     pos += (*riter) -> Thickness();
  }
 
  return pos;
}


void TargetGeometryManager::SetRadius(G4double rad) {

  layerRadius = rad;
  Notify();
}


void TargetGeometryManager::SetThickness(G4double thickn) {

  frontLayerThickness = thickn;
  Notify();
}


void TargetGeometryManager::SetMaterial(G4String mat) {

  frontLayerMaterial = mat;  
  Notify();
}

void TargetGeometryManager::SetMaxStepSize(G4double max) {

  maxStepSize = max;
  Notify();
}

