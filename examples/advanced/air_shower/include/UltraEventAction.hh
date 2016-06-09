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
// --------------------------------------------------------------
//                 GEANT 4 - ULTRA experiment example
// --------------------------------------------------------------
//
// Code developed by:
// B. Tome, M.C. Espirito-Santo, A. Trindade, P. Rodrigues 
//
//   **********************************************
//   *        UltraEventAction.hh
//   **********************************************
//
//    Ultra EventAction class. The UltraAnalysisManager class is used for histogram
//    filling if the G4ANALYSIS_USE environment variable is set.
//
#ifndef UltraEventAction_h
#define UltraEventAction_h 1

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  
#include "UltraRunAction.hh"

#include "G4UserEventAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

class UltraRunAction;
class UltraPrimaryGeneratorAction;
class G4HCofThisEvent;

class UltraEventAction : public G4UserEventAction
{
  public:
    UltraEventAction(UltraRunAction*);
   ~UltraEventAction();

  public:
    void BeginOfEventAction(const G4Event*);
    void EndOfEventAction(const G4Event*);
    

    G4int GetRunNumb(){return UltraRun->GetRunNumb() ;} ;
    G4int GetEvtNumb(){return evtNb ;} ;

  private:
    UltraRunAction* UltraRun;
    G4int    evtNb ;

  private: 

    G4int OpticalHitsCollID ; 
};
    
#endif
