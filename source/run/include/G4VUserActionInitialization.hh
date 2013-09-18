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
// $Id: G4VUserActionInitialization.hh 68830 2013-04-06 04:17:11Z asaim $
//

#ifndef G4VUserActionInitialization_h
#define G4VUserActionInitialization_h 1

// class description:
//
//  This is the abstract base class for instantiating all the user action classes.
// It has a pure virtual method Build() which is invoked by G4RunManager for
// sequential execution and G4WorkerRunManager for multi-threaded execution.
// The additional virtual method BuildForMaster() will be invoked from G4MTRunManager
// for multi-threaded execution.
//
//  Note that these virtual methods are const. It means the user may construct
// user action objects, but should not store the pointers of these objects as
// data members of the derived class.
//
// Note for multi-threaded mode:
//  The only action class the user may set to G4MTRunManager is a run action. It is
// then used at the beginning and the end of a run. It may be the same class or a
// dedicated class different from the run action instantiated for G4WorkerRunManager.
//

class G4VUserPrimaryGeneratorAction;
class G4UserRunAction;
class G4UserEventAction;
class G4UserStackingAction;
class G4UserTrackingAction;
class G4UserSteppingAction;
class G4VSteppingVerbose;

class G4VUserActionInitialization
{
  public:
    G4VUserActionInitialization();
    virtual ~G4VUserActionInitialization();

  public: // with description
    virtual void Build() const = 0;
    // Virtual method to be implemented by the user to instantiate user action
    // class objects. 
    virtual void BuildForMaster() const;
    // Virtual method to be implemented by the user to instantiate user run action
    // class object to be used by G4MTRunManager. This method is not invoked in
    // the sequential mode. The user should not use this method to instantiate
    // user action classes rather than user run action.
    virtual G4VSteppingVerbose* InitializeSteppingVerbose() const;
    // Virtual method to be implemented by the user if (s)he has a concrete
    // SteppingVerbose class to be used by the worker thread. In this case
    // (s)he should instantiate her/his SteppingVerbose in the concrete
    // implementation of this method and return its pointer. If this method is
    // not implemented, the default G4SteppingVerbose will be used. Please note
    // that this method affects only for the worker thread.

  protected: // with description
    void SetUserAction(G4VUserPrimaryGeneratorAction*) const;
    void SetUserAction(G4UserRunAction*) const;
    void SetUserAction(G4UserEventAction*) const;
    void SetUserAction(G4UserStackingAction*) const;
    void SetUserAction(G4UserTrackingAction*) const;
    void SetUserAction(G4UserSteppingAction*) const;
    // These methods should be used to define user's action classes.

};

#endif

