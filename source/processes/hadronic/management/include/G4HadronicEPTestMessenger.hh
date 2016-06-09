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
#ifndef G4HadronicEPTestMessenger_h
#define G4HadronicEPTestMessenger_h 1

// Class description:
// Messenger class to enable control of energy/momentum testing for
// hadronic processes and models.  

// Class Description - End

//
// Author: Dennis Wright (SLAC)
// Date:   1 April 2010
//

#include "G4UImessenger.hh"
#include "G4HadronicProcessStore.hh"

class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADouble;


class G4HadronicEPTestMessenger: public G4UImessenger
{
 public: //with description
   G4HadronicEPTestMessenger(G4HadronicProcessStore* theProcessStore);

   ~G4HadronicEPTestMessenger();

   void SetNewValue (G4UIcommand *command, G4String newValues);

 private:
   G4HadronicProcessStore* theProcessStore;
  
   G4UIdirectory* heptstDirectory;
   G4UIcmdWithAnInteger* reportLvlCmd;
   G4UIcmdWithADouble* procRelLvlCmd;
   G4UIcmdWithADoubleAndUnit* procAbsLvlCmd;
};

#endif

