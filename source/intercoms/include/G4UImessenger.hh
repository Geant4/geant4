// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UImessenger.hh,v 1.3 1999-11-11 15:36:05 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef G4UImessenger_h
#define G4UImessenger_h 1

#include "globals.hh"
#include "G4ios.hh"

class G4UIcommand;

// class description:
//
//  This class is the base class which represents a messenger which maintains
// the commands. The user who wants to define some commands must create his/her
// own concrete class derived from this class. The user's concrete messenger
// must have a responsibility of creating and deleting commands. Also, it must
// take care the delivering of the command to the destination class and replying
// the current value(s) of the parameter(s).
//

class G4UImessenger 
{

  public: // with description
      G4UImessenger();
      // Constructor. In the implementation of the concrete messenger, all commands
      // related to the messenger must be constructed.
      virtual ~G4UImessenger();
      // Destructor. In the implementation of the concrete messenger, all commands
      // defined in the constructor must be deleted.
      virtual G4String GetCurrentValue(G4UIcommand * command);
      // The concrete implementation of this method gets the current value(s) of the
      // parameter(s) of the given command from the destination class, converts the
      // value(s) to a string, and returns the string. Conversion could be done by
      // the ConvertToString() method of corresponding G4UIcmdXXX classes if the
      // the command is an object of these G4UIcmdXXX classes.
      virtual void SetNewValue(G4UIcommand * command,G4String newValue);
      // The concrete implementation of this method converts the string "newValue"
      // to value(s) of type(s) of the parameter(s). Convert methods corresponding
      // to the type of the command can be used if the command is an object of
      // G4UIcmdXXX classes.
  public:
      // For G4RWTPtrOrderedVector...
      G4bool operator == (const G4UImessenger& messenger) const;

  protected:
      G4String ItoS(G4int i);
      G4String DtoS(G4double a);
      G4String BtoS(G4bool b);
      G4int    StoI(G4String s);
      G4double StoD(G4String s);
      G4bool   StoB(G4String s);

  protected:
      void AddUIcommand(G4UIcommand * newCommand);
};

#endif

