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
// $Id: F03EventAction.hh,v 1.3 2001-11-07 16:36:33 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef F03EventAction_h
#define F03EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class F03RunAction;
class F03EventActionMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class F03EventAction : public G4UserEventAction
{
  public:
    F03EventAction(F03RunAction* F03RA);
   ~F03EventAction();

  public:
    void BeginOfEventAction(const G4Event*);
    void   EndOfEventAction(const G4Event*);
    G4int GetEventno();
    void setEventVerbose(G4int level);
    
    void SetDrawFlag(G4String val)  {drawFlag = val;};
    void SetPrintModulo(G4int val)  {printModulo = val;};
        
  private:
    G4int    calorimeterCollID;
    F03EventActionMessenger*  eventMessenger;
    F03RunAction* runaction;
    G4int verboselevel;
   
    G4String drawFlag;
    G4int    printModulo;             
};

#endif
