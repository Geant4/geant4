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
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4EnergyLossMessengerSTD
//
// Author:        
// 
// Creation date: 
//
// Modifications: 03.01.2002 V.Ivanchenko update to new design
//
//
// Class Description: 
//  This is a messenger class to interface to exchange information
//  between G4VEnergyLoss and UI.
//        
//  /process/eLoss/   directory
//
//   Commands :
//
//    rndmStep *     Randomize the proposed step by eLoss (false/true)
//    fluct *        Switch true/false the energy loss fluctuations
//    subsec *       Switch true/false the subcutoff generation
//    minsubsec *    Set the min. cut for subcutoff delta in range
//    StepFunction * Set the energy loss step limitation parameters
// 

// -------------------------------------------------------------------
//

#ifndef G4EnergyLossMessengerSTD_h
#define G4EnergyLossMessengerSTD_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class G4UIdirectory;
class G4UIcommand;
class G4UIcmdWithABool;
class G4UIcmdWithADoubleAndUnit;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4EnergyLossMessengerSTD: public G4UImessenger
{
public:   // with description
  
  G4EnergyLossMessengerSTD();
 ~G4EnergyLossMessengerSTD();
    
  void SetNewValue(G4UIcommand*, G4String);
    
private:

  G4EnergyLossMessengerSTD & operator=(const  G4EnergyLossMessengerSTD &right);
  G4EnergyLossMessengerSTD(const  G4EnergyLossMessengerSTD&);

  G4UIdirectory*             eLossDirectory;     
  G4UIcmdWithABool*          RndmStepCmd;
  G4UIcmdWithABool*          EnlossFlucCmd;
  G4UIcmdWithABool*          SubSecCmd;
  G4UIcmdWithABool*          IntegCmd;
  G4UIcmdWithADoubleAndUnit* MinSubSecCmd;
  G4UIcmdWithADoubleAndUnit* MinEnCmd;
  G4UIcmdWithADoubleAndUnit* MaxEnCmd;
  G4UIcommand*               StepFuncCmd;

};

#endif
