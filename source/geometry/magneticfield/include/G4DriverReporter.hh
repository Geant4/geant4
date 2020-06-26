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
// G4DriverReporter
//
// Class description:
//
// Auxiliary class to print information from integration drivers
//   Can be used by different types of drivers.

// Authors: J.Apostolakis        -   January/March 2020
// -------------------------------------------------------------------

#ifndef G4DRIVERREPORTER_HH
#define G4DRIVERREPORTER_HH

#include "G4FieldTrack.hh"

class G4DriverReporter 
{
  public:
    static void PrintStatus(const G4double* StartArr,
                            G4double        xstart,
                            const G4double* CurrentArr,
                            G4double        xcurrent,
                            G4double        requestStep,
                            unsigned int    subStepNo,
                            unsigned int    noIntegrationVariables);    
   
    static void PrintStatus(const G4FieldTrack& StartFT,
                            const G4FieldTrack& CurrentFT,
                            G4double            requestStep,
                            unsigned int        subStepNo);
   
    static void PrintStat_Aux(const G4FieldTrack& aFieldTrack,
                              G4double requestStep,
                              G4double actualStep,
                              G4int subStepNo,
                              G4double subStepSize,
                              G4double dotVelocities);
  
  private: 
    // G4int          fVerboseLevel;      // Verbose output for debugging
    // unsigned int   fNoIntegrationVariables);    
};
#endif
