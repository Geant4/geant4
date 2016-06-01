// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UserEventAction.hh,v 2.1 1998/07/12 02:53:56 urbi Exp $
// GEANT4 tag $Name: geant4-00 $
//
//
//

#ifndef G4UserEventAction_h
#define G4UserEventAction_h 1

class G4EventManager;

class G4UserEventAction 
{
  public:
      virtual ~G4UserEventAction() 
      {;}
      inline void SetEventManager(G4EventManager* value)
      { fpEventManager = value; };
      virtual void BeginOfEventAction();
      virtual void EndOfEventAction();
  protected:
      G4EventManager* fpEventManager;
};

#endif

