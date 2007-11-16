// $Id: EventAction.cc,v 1.1 2007-11-16 14:29:33 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   EventAction.cc
//
//                                         2007 Q
// ====================================================================
#include "EventAction.hh"
#include "G4Event.hh"
#include "Analysis.hh"

// ====================================================================
//
// class description
//
// ====================================================================

//////////////////////////
EventAction::EventAction()
//////////////////////////
{
}


///////////////////////////
EventAction::~EventAction()
///////////////////////////
{
}


////////////////////////////////////////////////////
void EventAction::BeginOfEventAction(const G4Event*)
////////////////////////////////////////////////////
{
}


//////////////////////////////////////////////////
void EventAction::EndOfEventAction(const G4Event*)
//////////////////////////////////////////////////
{
  Analysis* myana= Analysis::GetAnalysis();
  myana-> ClearIncidentFlag();
}

