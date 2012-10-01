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
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:		exrdm02RunActionMessenger.hh
//
// Author:		F Lei
// Organisation:	DERA UK
// Customer:		ESA/ESTEC, NOORDWIJK
// Contract:		12115/JG/NL Work Order No. 3
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// DESCRIPTION
// -----------

#ifndef exrdm02RunActionMessenger_h
#define exrdm02RunActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class exrdm02RunAction;
class G4UIcmdWithAString;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class exrdm02RunActionMessenger: public G4UImessenger
{
  public:
    exrdm02RunActionMessenger(exrdm02RunAction*);
   ~exrdm02RunActionMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    exrdm02RunAction*   exrdm02Run;   
    G4UIcmdWithAString* FileCmd;
};

#endif






