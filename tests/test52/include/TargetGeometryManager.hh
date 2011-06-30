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
#ifndef TARGETGEOMETRYMANAGER_HH
#define TARGETGEOMETRYMANAGER_HH

#include <vector>
#include "globals.hh"

class TargetComponent;
typedef std::vector<TargetComponent*> Components;


class TargetGeometryManager {

 public:
   TargetGeometryManager();
   ~TargetGeometryManager();

   void Attach(TargetComponent*);
   void Notify();

   TargetComponent* GetFrontLayer();
   G4double GetRadius() { return layerRadius; }
   G4double GetThickness(TargetComponent* comp);
   G4String GetMaterial(TargetComponent* comp);
   G4double GetPosition(TargetComponent* comp);
   G4double GetMaxStepSize() { return maxStepSize; }

   void SetRadius(G4double rad);
   void SetMaxStepSize(G4double max);
   void SetThickness(G4double thickn);
   void SetMaterial(G4String mat);

 private:
   Components* targetComponents;
    
   G4double layerRadius;
   G4double frontLayerThickness;
   G4String frontLayerMaterial;
   G4double maxStepSize;
};

#endif // TARGETGEOMETRYMANAGER_HH
