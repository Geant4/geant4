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
// $Id: G4MaterialScanner.hh 66892 2013-01-17 10:57:59Z gunter $
//
//


#ifndef G4MaterialScanner_H
#define G4MaterialScanner_H 1

// class description:
//
// G4MaterialScanner
//

#include "globals.hh"
#include "G4ThreeVector.hh"

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
  public: // with description
    G4MaterialScanner();

  public:
    ~G4MaterialScanner();

  public: // with description
    void Scan();
    // The main entry point which triggers ray tracing.
    // This method is available only if Geant4 is at Idle state.

  private:
    void DoScan();
    // Event loop
    void StoreUserActions();
    void RestoreUserActions();
    // Store and restore user action classes if defined

  private:
    G4RayShooter * theRayShooter;
    G4MatScanMessenger * theMessenger;

    G4EventManager * theEventManager;

    G4UserEventAction * theUserEventAction;
    G4UserStackingAction * theUserStackingAction;
    G4UserTrackingAction * theUserTrackingAction;
    G4UserSteppingAction * theUserSteppingAction;

    G4UserEventAction * theMatScannerEventAction;
    G4UserStackingAction * theMatScannerStackingAction;
    G4UserTrackingAction * theMatScannerTrackingAction;
    G4MSSteppingAction * theMatScannerSteppingAction;

    G4ThreeVector eyePosition;
    G4int nTheta;
    G4double thetaMin;
    G4double thetaSpan;
    G4int nPhi;
    G4double phiMin;
    G4double phiSpan;

    G4ThreeVector eyeDirection;

    G4bool regionSensitive;
    G4String regionName;
    G4Region* theRegion;

  public:
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
    inline void SetRegionSensitive(G4bool val=true) { regionSensitive = val; }
    inline G4bool GetRegionSensitive() const { return regionSensitive; }
    G4bool SetRegionName(const G4String& val);
    inline G4String GetRegionName() const { return regionName; }

};

#endif
