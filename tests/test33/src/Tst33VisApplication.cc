#include "Tst33VisApplication.hh"
#include "G4VisManager.hh"
#include "Tst33VisManager.hh"
#include "Tst33VisEventAction.hh"
#include "Tst33VisRunAction.hh"

#include "G4RunManager.hh"

#include "G4UImanager.hh"


Tst33VisApplication::  Tst33VisApplication() 
{
  fVisManager.Initialize();
}

Tst33VisApplication::~Tst33VisApplication(){
}

Tst33VEventAction *Tst33VisApplication::CreateEventAction() {
  return new Tst33VisEventAction;
}



G4UserRunAction *Tst33VisApplication::CreateRunAction(){
  return new Tst33VisRunAction;
}


