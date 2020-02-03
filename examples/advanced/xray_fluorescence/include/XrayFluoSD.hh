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
// Author: Elena Guardincerri (Elena.Guardincerri@ge.infn.it)
//
// History:
// -----------
//  28 Nov 2001  Elena Guardincerri   Created
//
// -------------------------------------------------------------------

#ifndef XrayFluoSD_hh
#define XrayFluoSD_hh 1

#include "G4VSensitiveDetector.hh"
#include "globals.hh"
//#include "XrayFluoDetectorConstruction.hh"
//#include "XrayFluoPlaneDetectorConstruction.hh"
//#include "XrayFluoMercuryDetectorConstruction.hh"
#include "XrayFluoSensorHit.hh"

class XrayFluoDetectorConstruction;
class XrayFluoPlaneDetectorConstruction;
class XrayFluoMercuryDetectorConstruction;
class G4HCofThisEvent;
class G4Step;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class XrayFluoSD : public G4VSensitiveDetector
{
public:
  
  XrayFluoSD(G4String, XrayFluoDetectorConstruction* );
  XrayFluoSD(G4String, XrayFluoPlaneDetectorConstruction* );
  XrayFluoSD(G4String, XrayFluoMercuryDetectorConstruction* );
  ~XrayFluoSD();
  
  void Initialize(G4HCofThisEvent*);
  void EndOfEvent(G4HCofThisEvent*);
  void clear();
  void DrawAll();
  void PrintAll();
  
protected:
  
 G4bool ProcessHits(G4Step*,G4TouchableHistory*);
private:
  
  //hit collection for the detector
  XrayFluoSensorHitsCollection*  HPGeCollection;   
  XrayFluoDetectorConstruction* Detector;
  XrayFluoPlaneDetectorConstruction* planeDetector;
  XrayFluoMercuryDetectorConstruction* mercuryDetector;

  G4int*                   HitHPGeID;
};

#endif









