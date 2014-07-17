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
// $Id: Tst69EventAction.hh 66512 2012-12-19 10:26:52Z gcosmo $
//
/////////////////////////////////////////////////////////////////////////
//
// Tst69EventAction
//
// Created: 17.07.2014 D. Mancusi
//
// Modified:
// 17.07.2014 Adapted from test46 (D. Mancusi)
//
////////////////////////////////////////////////////////////////////////
// 

#ifndef Tst69EventAction_h
#define Tst69EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include <ctime>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4Event;

class Tst69EventAction : public G4UserEventAction
{
  public:

    Tst69EventAction();
    virtual ~Tst69EventAction();

    void BeginOfEventAction(const G4Event*);
    void   EndOfEventAction(const G4Event*);

    G4int GetEventPrintModulo() const { return eventPrintModulo; }
    void SetEventPrintModulo(const G4int m) { eventPrintModulo = m; }

  private:
    clock_t actionStart;
    clock_t currentEventStart;
    clock_t lastEventStart;
    clock_t lastEventPrint;
    G4int lastEventIDPrint;
    G4int eventPrintModulo;
    G4int eventID;
    static const G4int nSecondsPerEventMax = 60;
    static const G4int nSecondsReport = 180;
};

#endif


