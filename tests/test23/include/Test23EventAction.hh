//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: Test23EventAction.hh,v 1.2 2006-06-29 21:49:31 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//      ---------------- Test23EventAction header ----------------
//                 by Mikhail Kossov, December 2003.
//  Test23EventAction class of the CHIPS Test of G4QCaptureAtRest process in GEANT4
// -----------------------------------------------------------------------------------
// It informs about the current event number 
// -----------------------------------------------------------------------------------

#ifndef Test23EventAction_h
#define Test23EventAction_h 1

#include "G4UserEventAction.hh"

#include "G4Event.hh"
#include "globals.hh"

class Test23EventActionMessenger;

class Test23EventAction : public G4UserEventAction
{
public: // Without description

    Test23EventAction();
   ~Test23EventAction();

    void BeginOfEventAction(const G4Event*);
    void   EndOfEventAction(const G4Event*);

    void SetDrawFlag(G4String val)  {drawFlag = val;};
    void SetPrintModulo(G4int val)  {printModulo = val;};

  private:

    G4int    nEvt;
    G4int    printModulo;
    G4int    verbose;
    G4String drawFlag;

    Test23EventActionMessenger*  eventMessenger;
};

#endif
