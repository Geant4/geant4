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
//   *        UltraRunAction.hh
//   **********************************************
//
//    RunAction class for Ultra; it has also a Messenger Class
//
#ifndef UltraRunAction_h
#define UltraRunAction_h 1

#include "globals.hh"
#include "G4UserRunAction.hh"

// forward declarations
class G4Run;
class UltraAnalysisManager;
class UltraRunActionMessenger;

class UltraRunAction : public G4UserRunAction
{
  public:

    UltraRunAction();
    ~UltraRunAction();

  public:

    void BeginOfRunAction(const G4Run* aRun);
    void EndOfRunAction(const G4Run* aRun);

    G4int GetRunNumb(){return runID;} ;

    void MySetRunID(G4int);

  private:
    UltraRunActionMessenger*     theRunActMessenger ;

    G4int saveRndm;
    G4int runID;
   
};

#endif 
