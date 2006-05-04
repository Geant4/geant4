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
// $Id: G4MaterialScanner.hh,v 1.1 2006-05-04 19:42:46 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
