// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UserEventAction.hh,v 1.4 1999-12-15 14:49:39 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//

#ifndef G4UserEventAction_h
#define G4UserEventAction_h 1

class G4EventManager;
class G4Event;

// class description:
//
//  This is the base class of one of the user's optional action classes.
// The two methods BeginOfEventAction() and EndOfEventAction() are invoked
// at the beginning and the end of one event processing. These methods are
// invoked by G4EventManager.
//  Be aware that BeginOfEventAction() is invoked when a G4Event object is
// sent to G4EventManager. Thus the primary vertexes/particles have already
// been made by the primary generator. In case the user wants to do something
// before generating primaries (i.e., store random number status), do it in
// the G4VUserPrimaryGeneratorAction concrete class.
//

class G4UserEventAction 
{
  public:
      G4UserEventAction() {;}
      virtual ~G4UserEventAction() {;}
      inline void SetEventManager(G4EventManager* value)
      { fpEventManager = value; }
  public: // with description
      virtual void BeginOfEventAction(const G4Event* anEvent);
      virtual void EndOfEventAction(const G4Event* anEvent);
      // Two virtual method the user can override.
  protected:
      G4EventManager* fpEventManager;
};

#endif

