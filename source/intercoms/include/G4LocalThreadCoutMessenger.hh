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
// $Id: G4LocalThreadCoutMessenger.hh 66241 2012-12-13 18:34:42Z gunter $
//
// class description
//
// This class is the messenger for handling cout/cerr of local thread
//

#ifndef G4LocalThreadCoutMessenger_h
#define G4LocalThreadCoutMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class G4UIdirectory;
class G4UIcommand;
class G4UIcmdWithABool;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4LocalThreadCoutMessenger: public G4UImessenger
{
  public:
    G4LocalThreadCoutMessenger();
   ~G4LocalThreadCoutMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:    
    G4UIdirectory*             coutDir;
    G4UIcommand*               coutFileNameCmd;
    G4UIcommand*               cerrFileNameCmd;
    G4UIcmdWithABool*      bufferCoutCmd;
    G4UIcmdWithAString*    prefixCmd;
    G4UIcmdWithAnInteger*  ignoreCmd;
    G4UIcmdWithABool*      ignoreInitCmd;
};

#endif

