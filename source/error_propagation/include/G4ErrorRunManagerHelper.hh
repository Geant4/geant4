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
// $Id: G4ErrorRunManagerHelper.hh 108476 2018-02-15 14:21:25Z gcosmo $
//
//
// Class Description:
//
// This class helps G4ErrorPropagatorManager with the initialization of 
// geometry and physics, definition of actions and run termination
// (things done by G4RunManager in an standard GEANT4 run)
// It holds a pointer to G4RunManagerKernel

// History:
// - Created:  Pedro Arce, January 2005
// --------------------------------------------------------------------

#ifndef G4ErrorRunManagerHelper_hh
#define G4ErrorRunManagerHelper_hh

#include "tls.hh"

class G4RunManagerKernel;
class G4VUserDetectorConstruction;
class G4VPhysicalVolume;
class G4VUserPhysicsList;
class G4UserTrackingAction;
class G4UserSteppingAction;

class G4ErrorRunManagerHelper 
{

 public:  // with description

  G4ErrorRunManagerHelper();
  virtual ~G4ErrorRunManagerHelper();

  static G4ErrorRunManagerHelper* GetRunManagerKernel();
    // Static method which returns the singleton pointer of
    // G4ErrorRunManagerHelper 

  void SetUserInitialization(G4VUserDetectorConstruction* userInit);
    // Initialize geometry by passing a G4VUserDetectorConstruction
  void SetUserInitialization(G4VPhysicalVolume* userInit);
    // Initialize geometry by passing a G4VPhysicalVolume

  void SetUserInitialization(G4VUserPhysicsList* userInit);
    // Initializes physics

  void SetUserAction(G4UserTrackingAction* userAction);
    // Set the user tracking action
  void SetUserAction(G4UserSteppingAction* userAction);
    // Set the user stepping action

  void RunInitialization();
    // Invokes G4RunManagerKernel RunInitialization();

  void InitializeGeometry();
    // Initializes GEANT4 geometry
  void InitializePhysics();
    // Initializes GEANT4 physics

  void RunTermination();
    // Invokes G4RunManagerKernel RunTermination();

  G4VUserPhysicsList* GetUserPhysicsList() const
    { return theUserPhysicsList; }

 private:

  static G4ThreadLocal G4ErrorRunManagerHelper* fRunManagerKernel;
  
  G4VUserPhysicsList* theUserPhysicsList;
  G4VPhysicalVolume* theUserWorld;

  G4RunManagerKernel* theG4RunManagerKernel;
};

#endif
