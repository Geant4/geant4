// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UImessenger.hh,v 1.1 1999-01-07 16:09:23 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef G4UImessenger_h
#define G4UImessenger_h 1

#include "globals.hh"
#include "G4ios.hh"

class G4UIcommand;

class G4UImessenger 
{

  public:
      G4UImessenger();
      virtual ~G4UImessenger();
      virtual G4String GetCurrentValue(G4UIcommand * command);
      virtual void SetNewValue(G4UIcommand * command,G4String newValue);
      // For RWTPtrOrderedVector...
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

