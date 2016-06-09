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
// $Id: HadrontherapyPrimaryGeneratorMessenger.hh; May 2005
// ----------------------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// ----------------------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone(a)*, F. Di Rosa(a), S. Guatelli(b), G. Russo(a)
// 
// (a) Laboratori Nazionali del Sud 
//     of the National Institute for Nuclear Physics, Catania, Italy
// (b) National Institute for Nuclear Physics Section of Genova, genova, Italy
// 
// * cirrone@lns.infn.it
// ----------------------------------------------------------------------------

#ifndef HadrontherapyPrimaryGeneratorMessenger_h
#define HadrontherapyPrimaryGeneratorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class HadrontherapyPrimaryGeneratorAction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithADouble;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class HadrontherapyPrimaryGeneratorMessenger: public G4UImessenger
{
  public:
    HadrontherapyPrimaryGeneratorMessenger(HadrontherapyPrimaryGeneratorAction*);
   ~HadrontherapyPrimaryGeneratorMessenger();
    
  void SetNewValue(G4UIcommand*, G4String);


 private:
    HadrontherapyPrimaryGeneratorAction* HadrontherapyAction; 
    G4UIdirectory*                    beamParametersDir;
    G4UIdirectory*                    EnergyDir;
    G4UIdirectory*                    particlePositionDir;
    G4UIdirectory*                    MomentumDir;
    G4UIcmdWithADoubleAndUnit*        meanKineticEnergyCmd;    
    G4UIcmdWithADoubleAndUnit*        sigmaEnergyCmd;  
    G4UIcmdWithADoubleAndUnit*        XpositionCmd;   
    G4UIcmdWithADoubleAndUnit*        YpositionCmd; 
    G4UIcmdWithADoubleAndUnit*        ZpositionCmd; 
    G4UIcmdWithADoubleAndUnit*        sigmaYCmd; 
    G4UIcmdWithADoubleAndUnit*        sigmaZCmd; 
    G4UIcmdWithADouble*        sigmaMomentumYCmd; 
    G4UIcmdWithADouble*        sigmaMomentumZCmd; 


};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

