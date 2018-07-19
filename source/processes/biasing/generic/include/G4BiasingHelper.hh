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
// $Id: $
//
//---------------------------------------------------------------
//
// G4BiasingHelper
//
// Class Description:
//    A utility class to help configuring code for biasing.
//    
//
//---------------------------------------------------------------
//   Initial version                         Sep. 2013 M. Verderi


#ifndef G4BiasingHelper_h
#define G4BiasingHelper_h 1

#include "globals.hh"

class G4ProcessManager;
class G4ParallelGeometriesLimiterProcess;

#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4BiasingHelper
{
public:
  // -- Substitute, in the process manager passed, physics process "physicsProcessToBias"
  // -- with a version wrapped by a G4BiasingProcessInterface wrapped, to put this physics
  // -- process under biasing.
  // -- Only processes of process type 2 (EM), 3 (Optical), 4 (Had.) and 6 (Decay) will be
  // -- wrapped.
  // -- A name for this wrapped process can be given (otherwise a default one is provided:
  // -- e.g. for process "phot" this will be "biasWrapper(phot)").
  static G4bool    ActivatePhysicsBiasing(G4ProcessManager* pmanager, G4String physicsProcessToBias,
					  G4String wrappedName = "");
  // -- Insert, in the process manager passed, a G4BiasingProcessInterface process that
  // -- will deal with non-modifying physics biasing (splitting, killing). A name for
  // -- this process can be passed, otherwise the default name "biasWrapper(0)" is used.
  static void   ActivateNonPhysicsBiasing(G4ProcessManager* pmanager,
					  G4String nonPhysicsProcessName = "");
  
  // -- Add a G4ParallelGeometriesLimiterProcess instance to the given process manager
  // -- The pointer of the added process is returned.
  static G4ParallelGeometriesLimiterProcess* AddLimiterProcess(G4ProcessManager* pmanager,
							       const G4String& processName = "biasLimiter");
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
