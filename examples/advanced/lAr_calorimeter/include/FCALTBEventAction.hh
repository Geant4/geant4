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
// $Id: FCALTBEventAction.hh,v 1.10 2004/11/29 18:03:06 ribon Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef FCALTBEventAction_h
#define FCALTBEventAction_h 1

#include "G4UserEventAction.hh"
#include "FCALSteppingAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

class FCALRunAction;
class FCALTBEventActionMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class FCALTBEventAction : public G4UserEventAction
{
  public:
    FCALTBEventAction(FCALSteppingAction*);
    virtual ~FCALTBEventAction();

  public:
    virtual void   BeginOfEventAction(const G4Event*);
    virtual void   EndOfEventAction(const G4Event*);
    
    void SetDrawFlag   (G4String val)  {drawFlag = val;};
    void SetPrintModulo(G4int    val)  {printModulo = val;};

  private:
    G4int                       calorimeterCollID;                
    G4String                    drawFlag;
    G4int                       printModulo;   
     
    FCALSteppingAction* StepAction;
    FCALTBEventActionMessenger*  eventMessenger;
    FCALRunAction* runManager;

  private:
  G4double NTracksOutOfWorld, NSecondaries, Init1, Init2, Init3;
  };

#endif    
