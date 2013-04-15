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
// $Id$
//
// Author: Ivana Hrivnacova, 10/04/2013  (ivana@ipno.in2p3.fr)
//
// The new user application base class which provides
// factory methods for instantiation of user application classes

#ifndef G4VUserApplication_hh
#define G4VUserApplication_hh

class G4VUserDetectorConstruction;
class G4VUserPrimaryGeneratorAction;
class G4VUserPhysicsList;
class G4UserRunAction;
class G4UserEventAction;
class G4UserStackingAction;
class G4UserTrackingAction;
class G4UserSteppingAction;

class G4VUserApplication {

  public:
    G4VUserApplication() {}
    virtual ~G4VUserApplication() {}
    
    // Create new instance
    virtual G4VUserApplication* CreateInstance() = 0;
    
    // User mandatory classes
    virtual G4VUserDetectorConstruction*   CreateDetectorConstruction() = 0;
    virtual G4VUserPhysicsList*            CreatePhysicsList() = 0;
    virtual G4VUserPrimaryGeneratorAction* CreatePrimaryGeneratorAction() = 0;
    
    // User action classes (optional)
    virtual G4UserRunAction*       CreateRunAction()      { return 0; }
    virtual G4UserEventAction*     CreateEventAction()    { return 0; }
    virtual G4UserStackingAction*  CreateStackingAction() { return 0; }
    virtual G4UserTrackingAction*  CreateTrackingAction() { return 0; }
    virtual G4UserSteppingAction*  CreateSteppingAction() { return 0; }

    // Function which can be used for additional settings
    virtual void Initialize() {}
};

#endif //G4VUserApplication_hh

