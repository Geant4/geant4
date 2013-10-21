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

class UltraRunAction : public G4UserRunAction
{
  public:

  UltraRunAction();
  ~UltraRunAction();
  
public:
  
  void BeginOfRunAction(const G4Run* aRun);
  void EndOfRunAction(const G4Run* aRun);  

private:  
  G4int saveRndm;
  G4int luxury;
  G4int seed;
   
};

#endif 
