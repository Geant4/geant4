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
// $Id: BrachyRunMessenger.hh,v 1.2 2002-11-18 15:18:37 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
//
//
//    ********************************************
//    *                                          *
//    *      BrachyRunMessenger.hh               *
//    *                                          *
//    ********************************************
// This class permits to switch the energy of the gamma delivered from the 
//radionuclides (Iodium/Iridium)
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef BrachyRunMessenger_h
#define BrachyRunMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class BrachyRunAction;
class BrachyFactoryIr;
class BrachyRunAction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class BrachyRunMessenger: public G4UImessenger
{
  public:
    BrachyRunMessenger(BrachyRunAction* );
   ~BrachyRunMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
  
  private:
 
     BrachyRunAction*  pRun;
     G4UIcmdWithAString*        selDetCmd;
    
    G4UIdirectory*             detDir;
  G4UIdirectory*           mydetDir;
    
  
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

