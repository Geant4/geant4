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
// $Id: Em0EventAction.hh,v 1.4 2001-07-11 10:03:36 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Em0EventAction_h
#define Em0EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class Em0EventActionMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Em0EventAction : public G4UserEventAction
{
  public:
    Em0EventAction();
   ~Em0EventAction();

  public:
    void BeginOfEventAction(const G4Event*);
    void   EndOfEventAction(const G4Event*);
    
    void addEdep(G4double Edep)    {TotalEnergyDeposit += Edep;};      
    G4double GetEnergyDeposit()    {return TotalEnergyDeposit;};    
    void SetDrawFlag(G4String val) {drawFlag = val;};
    
  private: 
    G4double                  TotalEnergyDeposit;   
    G4String                  drawFlag;             
    Em0EventActionMessenger*  eventMessenger;
};

#endif

    
