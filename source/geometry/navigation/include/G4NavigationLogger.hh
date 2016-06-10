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
// $Id: G4NavigationLogger.hh 93289 2015-10-15 10:01:15Z gcosmo $
//
// 
// class G4NavigationLogger
//
// Class description:
//
// Simple utility class for use by navigation systems
// for verbosity and check-mode.

// History:
// - Created. Gabriele Cosmo, November 2010
// --------------------------------------------------------------------
#ifndef G4NAVIGATIONLOGGER_HH
#define G4NAVIGATIONLOGGER_HH

#include "G4NavigationHistory.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4VSolid.hh"
#include "G4ThreeVector.hh"

class G4NavigationLogger
{
  public:  // with description

    G4NavigationLogger(const G4String& id);
   ~G4NavigationLogger();

    void PreComputeStepLog  (const G4VPhysicalVolume* motherPhysical,
                                   G4double motherSafety,
                             const G4ThreeVector& localPoint) const;
      // Report about first check - mother safety
   
    void AlongComputeStepLog(const G4VSolid* sampleSolid,
                             const G4ThreeVector& samplePoint,
                             const G4ThreeVector& sampleDirection,
                             const G4ThreeVector& localDirection,
                                   G4double sampleSafety,
                                   G4double sampleStep) const;
     // Report about a candidate daughter 
   
   void CheckDaughterEntryPoint(const G4VSolid* sampleSolid,
                                const G4ThreeVector& samplePoint,
                                const G4ThreeVector& sampleDirection,
                                const G4VSolid* motherSolid,
                                const G4ThreeVector& localPoint,
                                const G4ThreeVector& localDirection,
                                G4double motherStep,
                                G4double sampleStep) const;
     // Check suspicious distance to a candidate daughter
   
    void PostComputeStepLog (const G4VSolid* motherSolid,
                             const G4ThreeVector& localPoint,
                             const G4ThreeVector& localDirection,
                                   G4double motherStep,
                                   G4double motherSafety) const;
     // Report exit distance from mother 
   
    void ComputeSafetyLog   (const G4VSolid* solid,
                             const G4ThreeVector& point,
                                   G4double safety,
                                   G4bool isMotherVolume,    //  For labeling
                                   G4int banner= -1) const;
     // Report about safety computation (daughter?)
   
    void PrintDaughterLog   (const G4VSolid* sampleSolid,
                             const G4ThreeVector& samplePoint,
                                   G4double sampleSafety,
                                   G4bool   onlySafety,
                             const G4ThreeVector& sampleDirection,
                                   G4double sampleStep ) const;
     // Report about a new minimum distance to candidate daughter

    G4bool CheckAndReportBadNormal(const G4ThreeVector& unitNormal,
                                   const G4ThreeVector& localPoint,
                                   const G4ThreeVector& localDirection,
                                         G4double       step,
                                   const G4VSolid*      solid,                                 
                                   const char* msg ) const;
      // Report issue with normal from Solid    - for ComputeStep()

   G4bool CheckAndReportBadNormal(const G4ThreeVector& unitNormal,
                                  const G4ThreeVector& originalNormal,
                                  const G4RotationMatrix& rotationM,
                                  const char*          msg ) const;
      // Report issue with normal from Rotation - for ComputeStep()
   
    void ReportOutsideMother(const G4ThreeVector& localPoint,
                             const G4ThreeVector& localDirection,
                             const G4VPhysicalVolume* motherPV,
                                   G4double tDist = 30.0*CLHEP::cm ) const;
      // Report if point wrongly located outside mother volume

   void ReportVolumeAndIntersection( std::ostream& ostrm,
                                     const G4ThreeVector& localPoint,
                                     const G4ThreeVector& localDirection,
                                     const G4VPhysicalVolume* physical ) const;
      // Auxiliary method to report information about volume & position/direction 
      
  public:  // without description

    inline G4int GetVerboseLevel() const  { return fVerbose; }
    inline void  SetVerboseLevel(G4int level)  { fVerbose = level; }

    inline G4double GetMinTriggerDistance() const { return fMinTriggerDistance; }
    inline void     SetMinTriggerDistance(G4double d) { fMinTriggerDistance= d; }
    inline G4bool   GetReportSoftWarnings() const { return fReportSoftWarnings; }    
    inline void     SetReportSoftWarnings(G4bool b) { fReportSoftWarnings = b; } 
  private:

    G4String fId;                  // Navigation type
    G4int    fVerbose;             // Verbosity level
    G4double fMinTriggerDistance;  // Errors beyond this distance are fatal in ReportMother
    G4bool   fReportSoftWarnings;  // Flag> Whether to warn about small issues
};

#endif
