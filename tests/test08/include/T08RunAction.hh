// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: T08RunAction.hh,v 1.3 1999-04-22 22:09:58 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#ifndef T08RunAction_h
#define T08RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

class G4Run;

class T08RunAction : public G4UserRunAction
{
  public:
    T08RunAction();
    virtual ~T08RunAction();

  public:
    virtual void BeginOfRunAction(const G4Run* aRun);
    virtual void EndOfRunAction(const G4Run* aRun);

};

#endif

