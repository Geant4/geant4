//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: QEventAction.hh,v 1.2 2006-06-04 21:35:59 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   QEventAction.hh
//
//                                         2005 Q
// ====================================================================
#ifndef Q_EVENT_ACTION_H
#define Q_EVENT_ACTION_H

#include "G4UserEventAction.hh"

// ====================================================================
//
// class definition
//
// ====================================================================

class QEventAction : public G4UserEventAction {
public:
  QEventAction();
  ~QEventAction();

  virtual void BeginOfEventAction(const G4Event* aevent);
  virtual void EndOfEventAction(const G4Event* aevent);

};

#endif
