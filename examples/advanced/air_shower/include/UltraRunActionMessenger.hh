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
//   *        UltraRunActionMessenger.hh
//   **********************************************
//
//    Messenger Class for UltraRunAction
//    Allows to set the run ID
//
#ifndef UltraRunActionMessenger_h
#define UltraRunActionMessenger_h 1

class G4UIdirectory;
class UltraRunAction;
#include "G4UImessenger.hh"
#include "globals.hh"

class UltraRunActionMessenger: public G4UImessenger
{
public:
    UltraRunActionMessenger(UltraRunAction*);
    ~UltraRunActionMessenger();

  public:
    void SetNewValue(G4UIcommand * command,G4String newValues);

  private:
    UltraRunAction * theRunAction;
    
  private: //commands
    G4UIdirectory *             runDirectory;
    G4UIcommand *               runIDCmd ;
};

#endif


