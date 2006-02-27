// $Id: QEventAction.cc,v 1.1 2006-02-27 10:05:24 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   QEventAction.cc
//
//                                         2005 Q
// ====================================================================
#include "QEventAction.hh"
#include "globals.hh"

// ====================================================================
//
// class description
//
// ====================================================================

////////////////////////////
QEventAction::QEventAction()
////////////////////////////
{
}


/////////////////////////////
QEventAction::~QEventAction()
/////////////////////////////
{
}


////////////////////////////////////////////////////////////
void QEventAction::BeginOfEventAction(const G4Event* aevent)
////////////////////////////////////////////////////////////
{
  G4cout << "QEventAction::BofEA is called." << G4endl;
}


//////////////////////////////////////////////////////////
void QEventAction::EndOfEventAction(const G4Event* aevent)
//////////////////////////////////////////////////////////
{
  G4cout << "QEventAction::EofEA is called." << G4endl;
}

