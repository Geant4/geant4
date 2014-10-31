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
// -------------------------------------------------------------
//      GEANT 4 class 
//
//      ---------- TestStoppingPhysics -------
//                by Julia Yarba, 14 May 2009 
// 
//    Modified:
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef TestStoppingPhysics_h
#define TestStoppingPhysics_h 1

#include "globals.hh"

class G4VProcess;
class G4ProcessManager;
#if defined (USE_MUCAPTURE) 
class G4MuonMinusCapturePhysics;
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class TestStoppingPhysics
{
public:

  explicit TestStoppingPhysics(G4int verbose=0);

  ~TestStoppingPhysics();

  G4VProcess* GetProcess(const G4String&, const G4String&);

private:

#if defined (USE_MUCAPTURE) 
  G4MuonMinusCapturePhysics* theMuonMinusCaptureConstructor;
#endif

  G4VProcess*       theProcess;
  G4ProcessManager* theProcessMan;

  G4int verboseLevel;

};

#endif
