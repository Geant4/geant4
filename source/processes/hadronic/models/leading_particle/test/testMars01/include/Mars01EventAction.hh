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
//
// $Id: Mars01EventAction.hh,v 1.1 2001-12-13 14:58:42 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef Mars01EventAction_h
#define Mars01EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include "g4std/fstream"

class G4Event;

class Mars01EventAction : public G4UserEventAction
{
public:
  Mars01EventAction();
  virtual ~Mars01EventAction();
  
public:
  virtual void BeginOfEventAction(const G4Event*);
  virtual void EndOfEventAction(const G4Event*);
  
private:
  G4std::ofstream TallyFile;

};

#endif

    
