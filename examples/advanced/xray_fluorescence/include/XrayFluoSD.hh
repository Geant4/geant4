//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: XrayFluoSD.hh
// GEANT4 tag $Name: xray_fluo-V03-02-00
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









