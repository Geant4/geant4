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
// $Id: G4GMocrenMessenger.hh 68043 2013-03-13 14:27:49Z gcosmo $
//
//
// Created:  Mar. 31, 2009  Akinori Kimura  
//
// UI command definition for gMocren-file driver.
//
#ifndef G4GMOCRENMESSENGER_HH
#define G4GMOCRENMESSENGER_HH 1

#include "G4UImessenger.hh"
#include <vector>

class G4UIdirectory;
class G4UIcmdWithABool;
class G4UIcmdWithAString;
class G4UIcommand;
class G4UIcmdWithoutParameter;


class G4GMocrenMessenger : public G4UImessenger {
    
public:
  G4GMocrenMessenger();
  virtual ~G4GMocrenMessenger();

  virtual G4String GetCurrentValue(G4UIcommand * command);
  virtual void SetNewValue(G4UIcommand * command, G4String newValue);
        
  virtual G4String getEventNumberSuffix();
  virtual G4bool appendGeometry();
  virtual G4bool addPointAttributes();
  virtual G4bool useSolids();
  virtual G4bool writeInvisibles();
  virtual G4String getVolumeName();
  virtual std::vector<G4String> getHitNames();
  virtual G4String getScoringMeshName();
  virtual std::vector<G4String> getHitScorerNames();
  virtual void list();
  virtual void getNoVoxels(G4int & nx, G4int & ny, G4int & nz) const;
  virtual G4bool getDrawVolumeGrid() {return kDrawVolumeGrid;}

private:            
  G4UIdirectory* kgMocrenDirectory;
        
  G4String suffix;
  G4UIcmdWithAString* setEventNumberSuffixCommand;
        
  G4bool geometry;
  G4UIcmdWithABool* appendGeometryCommand;

  G4bool pointAttributes;
  G4UIcmdWithABool* addPointAttributesCommand;

  G4bool solids;
  G4UIcmdWithABool* useSolidsCommand;

  G4bool invisibles;

  G4String kgMocrenVolumeName;
  G4UIcmdWithAString* kSetgMocrenVolumeNameCommand;

  std::vector<G4String> kgMocrenHitNames;
  G4UIcmdWithAString* kAddgMocrenHitNameCommand;
  G4UIcmdWithoutParameter * kResetgMocrenHitNameCommand;

  G4String kgMocrenScoringMeshName;
  G4UIcmdWithAString * kSetgMocrenScoringMeshNameCommand;

  std::vector<G4String> kgMocrenHitScorerNames;
  G4UIcmdWithAString * kAddgMocrenHitScorerNameCommand;
  G4UIcmdWithoutParameter * kResetgMocrenHitScorerNameCommand;

  G4int kgMocrenNoVoxels[3];
  G4UIcommand * kSetgMocrenNoVoxelsCommand;

  G4UIcmdWithoutParameter * kListgMocrenCommand;

  G4bool kDrawVolumeGrid;
  G4UIcmdWithABool* kDrawVolumeGridCommand;
};

#endif
