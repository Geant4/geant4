// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UserRunAction.hh,v 1.4 1999-12-15 14:53:52 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef G4UserRunAction_h
#define G4UserRunAction_h 1

class G4Run;

// class description:
//
//  This is the base class of a user's action class which defines the
// user's action at the begining and the end of each run. The user can
// override the following two methods but the user should not change 
// any of the contents of G4Run object.
//    virtual void BeginOfRunAction(const G4Run* aRun);
//    virtual void EndOfRunAction(const G4Run* aRun);
//  The user's concrete class derived from this class must be set to
// G4RunManager via G4RunManager::SetUserAction() method.
//

class G4UserRunAction
{
  public:
    G4UserRunAction();
    virtual ~G4UserRunAction();

  public:
    virtual void BeginOfRunAction(const G4Run* aRun);
    virtual void EndOfRunAction(const G4Run* aRun);
};

#endif


