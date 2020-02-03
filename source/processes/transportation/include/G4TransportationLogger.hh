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
// 
// class G4TransportationLogger
//
// Class description:
//
// Simple utility class for use by navigation systems
// for verbosity and check-mode.

// History:
// - Created. John Apostolakis,  June 2018
// --------------------------------------------------------------------
#ifndef G4TRANSPORTATIONLOGGER_HH
#define G4TRANSPORTATIONLOGGER_HH

#include "globals.hh"

class G4Track;
class G4Step;

class G4TransportationLogger
{
  public:  // with description

    G4TransportationLogger(const G4String& className, G4int verbosity);   
    G4TransportationLogger(const char*     className, G4int verbosity);
   ~G4TransportationLogger();

    // Provide report in case of particle looping in field
    //
    void ReportLoopingTrack( const G4Track & track,
                             const G4Step  & stepInfo,
                             G4int           numTrials,                            
                             G4long          noCalls,
                             const char*     methodName) const;

    // Print the thresholds' values
    void ReportLooperThresholds( const char* className );
   
  public:  // without description

    void SetThresholds( G4double newEnWarn, G4double importantEnergy,
                        G4int newMaxTrials );
   
    G4int GetVerboseLevel() const  { return fVerbose; }
    void  SetVerboseLevel(G4int level)  { fVerbose = level; }

    // Get/Set limit parameters for use in reporting
    //
    G4double GetThresholdWarningEnergy() const { return fThldWarningEnergy; }
    G4double GetThresholdImportantEnergy() const { return fThldImportantEnergy; }
    G4double GetThresholdTrials() const { return fThldTrials; }
       
    void SetThresholdWarningEnergy( G4double val )  { fThldWarningEnergy= val; }
    void SetThresholdImportantEnergy( G4double val ) { fThldImportantEnergy= val; }
    void SetThresholdTrials(G4int maxNoTrials ) { fThldTrials = std::max( maxNoTrials, 1); }
  
  private:           

    G4String fClassName;       // Name of Transportation process (class name)
    G4int    fVerbose;         // Verbosity level

    // Parameters for transporation limits
    // Used only for reporting in this class
    //
    G4double fThldWarningEnergy;
    G4double fThldImportantEnergy;
    G4int fThldTrials;
};

#endif
