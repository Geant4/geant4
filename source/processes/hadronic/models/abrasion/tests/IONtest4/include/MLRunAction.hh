#ifndef MLRunAction_h
#define MLRunAction_h 1
////////////////////////////////////////////////////////////////////////////////
//
#include "G4UserRunAction.hh"
#include "globals.hh"

#include "G4Run.hh"
////////////////////////////////////////////////////////////////////////////////
//
class MLRunAction : public G4UserRunAction
{
public:
  MLRunAction ();
  ~MLRunAction ();
  
public:
  void BeginOfRunAction (const G4Run*);
  void EndOfRunAction (const G4Run*);

};
////////////////////////////////////////////////////////////////////////////////
#endif
