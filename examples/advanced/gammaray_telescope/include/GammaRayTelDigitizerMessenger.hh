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
// $Id: GammaRayTelDigitizerMessenger.hh,v 1.1 2001-11-29 09:34:17 flongo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file
//      CERN Geneva Switzerland
//
//
//      ------------ GammaRayTelDigitizerMessenger  ------
//           by F.Longo, G.Santin & R.Giannitrapani (29 nov 2001) 
//
// ************************************************************

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef GammaRayTelDigitizerMessenger_h
#define GammaRayTelDigitizerMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class GammaRayTelDigitizer;
class G4UIcmdWithADoubleAndUnit;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class GammaRayTelDigitizerMessenger: public G4UImessenger
{
public:
  GammaRayTelDigitizerMessenger(GammaRayTelDigitizer*);
  ~GammaRayTelDigitizerMessenger();
  
  void SetNewValue(G4UIcommand*, G4String);
  
private:
  GammaRayTelDigitizer* GammaRayTelAction; 
  G4UIcmdWithADoubleAndUnit*  ThresholdCmd;

};

#endif


