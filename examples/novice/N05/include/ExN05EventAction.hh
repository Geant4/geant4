// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN05EventAction.hh,v 1.1 1999-01-07 16:06:12 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef ExN05EventAction_h
#define ExN05EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class G4Event;

class ExN05EventAction : public G4UserEventAction
{
  public:
    ExN05EventAction();
    ~ExN05EventAction();

  public:
    void BeginOfEventAction();
    void EndOfEventAction();

  private:
    G4int calorimeterCollID;
    G4int hadCalorimeterCollID;
    G4bool drawFlag;

  public:
    inline void SetDrawFlag(G4bool val)
    { drawFlag = val; };
};

#endif

    
