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
#ifndef G4NeutronHPMessenger_h
#define G4NeutronHPMessenger_h

#include "globals.hh"
#include "G4UImessenger.hh"

class G4NeutronHPManager;
class G4UIdirectory;
class G4UIcmdWithAString;

class G4NeutronHPMessenger: public G4UImessenger
{
   public:
      G4NeutronHPMessenger( G4NeutronHPManager* );
     ~G4NeutronHPMessenger();

      void SetNewValue(G4UIcommand*, G4String);

   private:
      G4NeutronHPManager* manager;

      G4UIdirectory* NeutronHPDir;
      G4UIcmdWithAString* PhotoEvaCmd;
      G4UIcmdWithAString* SkipMissingCmd;
      G4UIcmdWithAString* NeglectDopplerCmd;
      G4UIcmdWithAString* DoNotAdjustFSCmd;
      G4UIcmdWithAString* ProduceFissionFragementCmd;
      //G4UIcmdWithAString* AllowHeavyElementCmd;
/*
 * #setenv G4NEUTRONHP_USE_ONLY_PHOTONEVAPORATION 1
 * #setenv G4NEUTRONHP_SKIP_MISSING_ISOTOPES 1 
 * #setenv G4NEUTRONHP_NEGLECT_DOPPLER 1
 * #setenv G4NEUTRONHP_DO_NOT_ADJUST_FINAL_STATE 1
 * #setenv G4NEUTRONHP_PRODUCE_FISSION_FRAGMENTS 1
 *
*/
    
};

#endif
