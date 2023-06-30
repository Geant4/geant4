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
// G4MaterialScanner
//
// Class description:
//
// Utility class for scanning of materials through ray-tracing
// in a detector setup.

// Author: M.Asai, May 2006
// --------------------------------------------------------------------
#ifndef G4MaterialScanner_hh
#define G4MaterialScanner_hh 1

#include "G4ThreeVector.hh"
#include "globals.hh"

class G4Event;
class G4EventManager;
class G4UserEventAction;
class G4UserStackingAction;
class G4UserTrackingAction;
class G4UserSteppingAction;
class G4MSSteppingAction;
class G4MatScanMessenger;
class G4RayShooter;
class G4Region;

class G4MaterialScanner
{
  public:
    G4MaterialScanner();
    ~G4MaterialScanner();

    // The main entry point which triggers ray tracing.
    // This method is available only if Geant4 is in Idle state.
    void Scan();

    inline void SetEyePosition(const G4ThreeVector& val) { eyePosition = val; }
    inline G4ThreeVector GetEyePosition() const { return eyePosition; }
    inline void SetNTheta(G4int val) { nTheta = val; }
    inline G4int GetNTheta() const { return nTheta; }
    inline void SetThetaMin(G4double val) { thetaMin = val; }
    inline G4double GetThetaMin() const { return thetaMin; }
    inline void SetThetaSpan(G4double val) { thetaSpan = val; }
    inline G4double GetThetaSpan() const { return thetaSpan; }
    inline void SetNPhi(G4int val) { nPhi = val; }
    inline G4int GetNPhi() const { return nPhi; }
    inline void SetPhiMin(G4double val) { phiMin = val; }
    inline G4double GetPhiMin() const { return phiMin; }
    inline void SetPhiSpan(G4double val) { phiSpan = val; }
    inline G4double GetPhiSpan() const { return phiSpan; }
    inline void SetRegionSensitive(G4bool val = true) { regionSensitive = val; }
    inline G4bool GetRegionSensitive() const { return regionSensitive; }
    G4bool SetRegionName(const G4String& val);
    inline const G4String& GetRegionName() const { return regionName; }

  private:
    void DoScan();
    // Event loop

    void StoreUserActions();
    void RestoreUserActions();
    // Store and restore user action classes if defined.

  private:
    G4RayShooter* theRayShooter = nullptr;
    G4MatScanMessenger* theMessenger = nullptr;

    G4EventManager* theEventManager = nullptr;

    G4UserEventAction* theUserEventAction = nullptr;
    G4UserStackingAction* theUserStackingAction = nullptr;
    G4UserTrackingAction* theUserTrackingAction = nullptr;
    G4UserSteppingAction* theUserSteppingAction = nullptr;

    G4UserEventAction* theMatScannerEventAction = nullptr;
    G4UserStackingAction* theMatScannerStackingAction = nullptr;
    G4UserTrackingAction* theMatScannerTrackingAction = nullptr;
    G4MSSteppingAction* theMatScannerSteppingAction = nullptr;

    G4ThreeVector eyePosition;
    G4int nTheta = 91;
    G4double thetaMin = 0.0;
    G4double thetaSpan = 0.0;
    G4int nPhi = 37;
    G4double phiMin = 0.0;
    G4double phiSpan = 0.0;

    G4ThreeVector eyeDirection;

    G4bool regionSensitive = false;
    G4String regionName = "notDefined";
    G4Region* theRegion = nullptr;
};

#endif
