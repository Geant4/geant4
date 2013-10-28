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
#ifndef TstTarget_H
#define TstTarget_H 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"       

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;


class TstTarget : public G4VUserDetectorConstruction {

public:

  TstTarget();
  ~TstTarget();
  
  G4VPhysicalVolume*         ConstructTarget();
  virtual G4VPhysicalVolume* Construct();
  void                       ResetGeom();
  
  void        SetDimentions( G4double x1, G4double x2, G4double x3 ) { fX1=x1; fX2=x2, fX3=x3; return; }
  void        SetX( G4double val ) { fX1=val; return; }
  void        SetY( G4double val ) { fX2=val; return; }
  void        SetZ( G4double val ) { fX3=val; return; }
  void        SetRMin( G4double val ) { fX1=val; return; }
  void        SetRMax( G4double val ) { fX2=val; return; }
  void        SetShape( G4String val ) { fShape=val; return; }
  
  G4Material* ResetMaterial( G4String mat );
  void        ResetMaterial( G4Material* mat ); 
  
  G4Material*        GetCurrentMaterial() const { return fMaterial; }
  G4VPhysicalVolume* GetTarget()          const { return fPhysTarget; }
  
private:

  G4double fX1;
  G4double fX2;
  G4double fX3;
  
  G4String fShape;
  
  G4String fMatName;
  G4Material* fMaterial;
  
  G4VPhysicalVolume* fWorld;
  G4VPhysicalVolume* fSubWorld;

  // the real stuff
  //
  G4LogicalVolume* fLogTarget;
  G4VPhysicalVolume* fPhysTarget;
  
};

#endif
