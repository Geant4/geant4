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
// Hadrontherapy advanced example for Geant4
// See more at: https://twiki.cern.ch/twiki/bin/view/Geant4/AdvancedExamplesHadrontherapy

#include "HadrontherapyPrimaryGeneratorMessenger.hh"
#include "HadrontherapyPrimaryGeneratorAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4SystemOfUnits.hh"

HadrontherapyPrimaryGeneratorMessenger::HadrontherapyPrimaryGeneratorMessenger(
                                                                               HadrontherapyPrimaryGeneratorAction* HadrontherapyGun)
:HadrontherapyAction(HadrontherapyGun)
<<<<<<< HEAD
{ 
  //
  // Definition of the interactive commands to modify the parameters of the
  // generation of primary particles
  // 
 beamParametersDir = new G4UIdirectory("/beam/");
 beamParametersDir -> SetGuidance("set parameters of beam");
 
 EnergyDir = new G4UIdirectory("/beam/energy/");  
 EnergyDir -> SetGuidance ("set energy of beam");  

 particlePositionDir = new G4UIdirectory("/beam/position/");  
 particlePositionDir -> SetGuidance ("set position of particle");  

 MomentumDir = new G4UIdirectory("/beam/momentum/");  
 MomentumDir -> SetGuidance ("set momentum of particle ");  
/*
 sigmaMomentumYCmd = new G4UIcmdWithADouble("/beam/momentum/sigmaY",this);
 sigmaMomentumYCmd -> SetGuidance("set sigma momentum y");
 sigmaMomentumYCmd -> SetParameterName("momentum",false);
 sigmaMomentumYCmd -> AvailableForStates(G4State_PreInit,G4State_Idle);   

 sigmaMomentumZCmd = new G4UIcmdWithADouble("/beam/momentum/sigmaZ",this);
 sigmaMomentumZCmd -> SetGuidance("set sigma momentum z");
 sigmaMomentumZCmd -> SetParameterName("momentum",false);
 sigmaMomentumZCmd -> AvailableForStates(G4State_PreInit,G4State_Idle);   

 meanKineticEnergyCmd = new G4UIcmdWithADoubleAndUnit("/beam/energy/meanEnergy",this);
 meanKineticEnergyCmd -> SetGuidance("set mean Kinetic energy");
 meanKineticEnergyCmd -> SetParameterName("Energy",false);
 meanKineticEnergyCmd -> SetDefaultUnit("MeV");
 meanKineticEnergyCmd -> SetUnitCandidates("eV keV MeV GeV TeV");
 meanKineticEnergyCmd -> AvailableForStates(G4State_PreInit,G4State_Idle);   
 
 sigmaEnergyCmd = new G4UIcmdWithADoubleAndUnit("/beam/energy/sigmaEnergy",this);
 sigmaEnergyCmd -> SetGuidance("set sigma energy");
 sigmaEnergyCmd -> SetParameterName("Energy",false);
 sigmaEnergyCmd -> SetDefaultUnit("keV");
 sigmaEnergyCmd -> SetUnitCandidates("eV keV MeV GeV TeV");
 sigmaEnergyCmd -> AvailableForStates(G4State_PreInit,G4State_Idle);   
 
 XpositionCmd = new G4UIcmdWithADoubleAndUnit("/beam/position/Xposition",this);
 XpositionCmd -> SetGuidance("set x coordinate of particle");
 XpositionCmd -> SetParameterName("position",false);
 XpositionCmd -> SetDefaultUnit("mm");
 XpositionCmd -> SetUnitCandidates("mm cm m");
 XpositionCmd -> AvailableForStates(G4State_PreInit,G4State_Idle);   

 YpositionCmd = new G4UIcmdWithADoubleAndUnit("/beam/position/Yposition",this);
 YpositionCmd -> SetGuidance("set y coordinate of particle");
 YpositionCmd -> SetParameterName("position",false);
 YpositionCmd -> SetDefaultUnit("mm");
 YpositionCmd -> SetUnitCandidates("mm cm m");
 YpositionCmd -> AvailableForStates(G4State_PreInit,G4State_Idle);   

 sigmaYCmd = new G4UIcmdWithADoubleAndUnit("/beam/position/Yposition/sigmaY",this);
 sigmaYCmd -> SetGuidance("set sigma y");
 sigmaYCmd -> SetParameterName("position",false);
 sigmaYCmd -> SetDefaultUnit("mm");
 sigmaYCmd -> SetUnitCandidates("mm cm m");
 sigmaYCmd -> AvailableForStates(G4State_PreInit,G4State_Idle);   

 ZpositionCmd = new G4UIcmdWithADoubleAndUnit("/beam/position/Zposition",this);
 ZpositionCmd -> SetGuidance("set z coordinate of particle");
 ZpositionCmd -> SetParameterName("position",false);
 ZpositionCmd -> SetDefaultUnit("mm");
 ZpositionCmd -> SetUnitCandidates("mm cm m");
 ZpositionCmd -> AvailableForStates(G4State_PreInit,G4State_Idle);   

 sigmaZCmd = new G4UIcmdWithADoubleAndUnit("/beam/position/Zposition/sigmaZ",this);
 sigmaZCmd -> SetGuidance("set sigma z");
 sigmaZCmd -> SetParameterName("position",false);
 sigmaZCmd -> SetDefaultUnit("mm");
 sigmaZCmd -> SetUnitCandidates("mm cm m");
 sigmaZCmd -> AvailableForStates(G4State_PreInit,G4State_Idle);   
*/
=======
{
    
    changeTheSource = new G4UIdirectory("/primaryParticleData/");
    changeTheSource -> SetGuidance("Command to change the source");
    
    NewSource = new G4UIcmdWithABool("/primaryParticleData/NewSource",this);
    NewSource -> SetParameterName("NewSource ", true);
    NewSource -> SetDefaultValue("false");
    NewSource -> SetGuidance("Set if you want read a .txt file"
                             "\n[usage]: /primaryParticleData/NewSource [true/false]");
    NewSource -> AvailableForStates(G4State_Idle, G4State_PreInit);
    
    
    calculatedPhaseSpaceFileIN=new G4UIcmdWithAString("/primaryParticleData/calculatedPhaseSpaceFileIN",this);
    calculatedPhaseSpaceFileIN->SetDefaultValue("");
    calculatedPhaseSpaceFileIN->SetGuidance("full path and file name of the phase space file to be used as particle generator");
    
    
    
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
}

HadrontherapyPrimaryGeneratorMessenger::~HadrontherapyPrimaryGeneratorMessenger()
{
<<<<<<< HEAD
  delete beamParametersDir;
  delete EnergyDir;
 // delete meanKineticEnergyCmd;  
 // delete sigmaEnergyCmd;
  delete particlePositionDir;
  delete MomentumDir;
/*  delete XpositionCmd; 
  delete YpositionCmd; 
  delete ZpositionCmd; 
  delete sigmaYCmd; 
  delete sigmaZCmd; 
  delete sigmaMomentumYCmd; 
  delete sigmaMomentumZCmd; */
}  
/*
void HadrontherapyPrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)                
{
  if ( command == meanKineticEnergyCmd )                                                                        
    { HadrontherapyAction -> SetmeanKineticEnergy(meanKineticEnergyCmd
						  -> GetNewDoubleValue(newValue));}
  if ( command == sigmaEnergyCmd )                                                                        
    { HadrontherapyAction -> SetsigmaEnergy(sigmaEnergyCmd
					    -> GetNewDoubleValue(newValue));}
  if ( command == XpositionCmd )                                                                        
    { HadrontherapyAction -> SetXposition(XpositionCmd
					  -> GetNewDoubleValue(newValue));}

  if ( command == YpositionCmd )                                                                        
    { HadrontherapyAction -> SetYposition(YpositionCmd
					  -> GetNewDoubleValue(newValue));}

  if ( command == ZpositionCmd )                                                                        
    { HadrontherapyAction -> SetZposition(ZpositionCmd
					  -> GetNewDoubleValue(newValue));}

  if ( command == sigmaYCmd )                                                                        
    { HadrontherapyAction -> SetsigmaY(sigmaYCmd
				       -> GetNewDoubleValue(newValue));}
=======
    
    delete NewSource;
    delete changeTheSource;
    delete calculatedPhaseSpaceFileIN;
    
    
}


void HadrontherapyPrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command, G4String NewValue)

{
    if (command == calculatedPhaseSpaceFileIN)
    {
        HadrontherapyAction->setCalculatedPhaseSpaceFileIN(NewValue);
    }
    
    if (command == NewSource)
    {
        HadrontherapyAction->setNewSource(NewSource->GetNewBoolValue(NewValue));
        
    }
    
}
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c

  if ( command == sigmaZCmd )                                                                        
    { HadrontherapyAction -> SetsigmaZ(sigmaZCmd
				       -> GetNewDoubleValue(newValue));}

  if ( command == sigmaMomentumYCmd )                                                                        
    { HadrontherapyAction -> SetsigmaMomentumY(sigmaMomentumYCmd
					       -> GetNewDoubleValue(newValue));}

  if ( command == sigmaMomentumZCmd )                                                                        
    { HadrontherapyAction -> SetsigmaMomentumZ(sigmaMomentumZCmd
					       -> GetNewDoubleValue(newValue));}
}                 */
