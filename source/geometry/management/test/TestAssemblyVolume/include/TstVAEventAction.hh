// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: TstVAEventAction.hh,v 1.3 2001-02-07 17:30:59 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// --------------------------------------------------------------

#ifndef TstVAEventAction_h
#define TstVAEventAction_h 1

#include "globals.hh"
#include "G4UserEventAction.hh"

class TstVAEventActionMessenger;

class TstVAEventAction : public G4UserEventAction
{
  public:
    TstVAEventAction();
    virtual ~TstVAEventAction();

  public:
    virtual void   BeginOfEventAction(const G4Event*);
    virtual void   EndOfEventAction(const G4Event*);
    
    void SetDrawFlag   (G4String val)  {drawFlag = val;};
    void SetPrintModulo(G4int    val)  {printModulo = val;};
    
  private:
    G4int                       calorimeterCollID;                
    G4String                    drawFlag;
    G4int                       printModulo;                         
    TstVAEventActionMessenger*  eventMessenger;
};

#endif

    
