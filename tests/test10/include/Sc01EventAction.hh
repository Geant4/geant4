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
// $Id: Sc01EventAction.hh,v 1.1 2004-01-27 11:19:21 grichine Exp $
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//      This class is a class derived from G4VUserEventAction
//      for constructing all particles and processes.
//
//	History
//        first version              09 Sept. 1998 by S.Magni
// ------------------------------------------------------------

#ifndef Sc01EventAction_h
#define Sc01EventAction_h 1

#include "G4UserEventAction.hh"

#include "globals.hh"

class Sc01EventActionMessenger;
class G4Event;

class Sc01EventAction : public G4UserEventAction
{
  public:
    Sc01EventAction();
    virtual ~Sc01EventAction();
    void   SetPrintModulo(G4int val)  { printModulo = val;};

  public:
    virtual void BeginOfEventAction(const G4Event*);
    virtual void EndOfEventAction(const G4Event*);

private:

  G4int    printModulo;
  Sc01EventActionMessenger*  eventMessenger;
};

#endif

    
