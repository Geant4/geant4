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
// $Id: Em5EventAction.hh,v 1.4 2001-10-16 11:56:27 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef Em5EventAction_h
#define Em5EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class Em5RunAction;
class Em5EventActionMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Em5EventAction : public G4UserEventAction
{
  public:
    Em5EventAction(Em5RunAction* Em5RA);
   ~Em5EventAction();

  public:
    void BeginOfEventAction(const G4Event*);
    void   EndOfEventAction(const G4Event*);
    G4int GetEventno();
    void setEventVerbose(G4int level);
    
    void CountStepsCharged() ;
    void CountStepsNeutral() ;
    void AddCharged() ;
    void AddNeutral() ;
    void AddE();
    void AddP();   
    void SetTr();
    void SetRef();
    
    void SetDrawFlag(G4String val)  {drawFlag = val;};
    void SetPrintModulo(G4int val)  {printModulo = val;};
        
  private:
    G4int    calorimeterCollID;
    Em5EventActionMessenger*  eventMessenger;
    Em5RunAction* runaction;
    G4int verboselevel;
    G4double nstep,nstepCharged,nstepNeutral;
    G4double Nch,Nne;
    G4double NE,NP;
    G4double Transmitted,Reflected ;
    
    G4String drawFlag;
    G4int    printModulo;             
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
