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
// $Id: G4AdjointPhysicsMessenger.cc,v 1.1 2009-11-19 22:41:18 ldesorgh Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//////////////////////////////////////////////////////////////
//      Class Name:	G4AdjointPhysicsMessenger
//	Author:       	L. Desorgher
// 	Organisation: 	SpaceIT GmbH
//	Contract:	ESA contract 21435/08/NL/AT
// 	Customer:     	ESA/ESTEC
//////////////////////////////////////////////////////////////
#include "G4AdjointPhysicsMessenger.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UnitsTable.hh"
#include "G4AdjointPhysicsList.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
G4AdjointPhysicsMessenger::G4AdjointPhysicsMessenger(
                                           G4AdjointPhysicsList* pPhysicsList)
:thePhysicsList(pPhysicsList)
{ 
  PhysicsDir = new G4UIdirectory("/adjoint_physics/");
  PhysicsDir->SetGuidance("Definition of physics processes to be used in the adjoint and forward simulation mode");
 
  
  
  


  //Physics
  //-------
  
 /* UseIonisationCmd = new G4UIcmdWithABool("/adjoint_physics/UseElectronIonisation",this);
  UseIonisationCmd->SetGuidance("If true (false) the ionisation is (not) considered for both adjoint and forward simulations");
  UseIonisationCmd->AvailableForStates(G4State_PreInit);
 */ 
  
  UsepIonisationCmd = new G4UIcmdWithABool("/adjoint_physics/UseProtonIonisation",this);
  UsepIonisationCmd->SetGuidance("If true (false) the proton ionisation is (not) considered for both adjoint and forward simulations");
  UsepIonisationCmd->AvailableForStates(G4State_PreInit);
  
  UseBremCmd = new G4UIcmdWithABool("/adjoint_physics/UseBremsstrahlung",this);
  UseBremCmd->SetGuidance("If true (false) the bremsstrahlung process is (not) considered for both adjoint and forward simulations");
  UseBremCmd->AvailableForStates(G4State_PreInit); 
  
  UseComptonCmd = new G4UIcmdWithABool("/adjoint_physics/UseCompton",this);
  UseComptonCmd->SetGuidance("If true (false) the Compton scattering is (not) considered for both adjoint and forward simulations");
  UseComptonCmd->AvailableForStates(G4State_PreInit);  
  
  UseMSCmd = new G4UIcmdWithABool("/adjoint_physics/UseMS",this);
  UseMSCmd->SetGuidance("If true (false) the continuous multiple scattering is (not) considered for both adjoint and forward simulations");
  UseMSCmd->AvailableForStates(G4State_PreInit); 
  
  UseEgainFluctuationCmd = new G4UIcmdWithABool("/adjoint_physics/UseEgainElossFluctuation",this);
  UseEgainFluctuationCmd->SetGuidance("If true (false) the fluctation for continuous energy gain and loss is (not) considered for adjoint and forward simulations");
  UseEgainFluctuationCmd->AvailableForStates(G4State_PreInit); 
  
  UsePEEffectCmd = new G4UIcmdWithABool("/adjoint_physics/UsePEEffect",this);
  UsePEEffectCmd->AvailableForStates(G4State_PreInit); 
  UsePEEffectCmd->SetGuidance("If true (false) the photo electric effect is (not) considered for the adjoint and forward simulations");
   
  UseGammaConversionCmd = new G4UIcmdWithABool("/adjoint_physics/UseGammaConversion",this);
  UseGammaConversionCmd->AvailableForStates(G4State_PreInit);
  UseGammaConversionCmd->SetGuidance("If true the gamma pair conversion is considered as weel as the e+ physics for the forward simulation");  
  
  SetEminAdjModelsCmd =  new G4UIcmdWithADoubleAndUnit("/adjoint_physics/SetEminForAdjointModels",this);
  SetEminAdjModelsCmd->SetGuidance("Set the minimum energy  of the adjoint models");
  SetEminAdjModelsCmd->SetParameterName("Emin",false);
  SetEminAdjModelsCmd->SetUnitCategory("Energy");
  SetEminAdjModelsCmd->AvailableForStates(G4State_PreInit);
  
  SetEmaxAdjModelsCmd =  new G4UIcmdWithADoubleAndUnit("/adjoint_physics/SetEmaxForAdjointModels",this);
  SetEmaxAdjModelsCmd->SetGuidance("Set the minimum energy  of the adjoint models.");
  SetEmaxAdjModelsCmd->SetParameterName("Emax",false);
  SetEmaxAdjModelsCmd->SetUnitCategory("Energy");
  SetEmaxAdjModelsCmd->AvailableForStates(G4State_PreInit);
/*
#ifdef TEST_MODE   
  SetCSBiasingFactorComptonCmd =  new G4UIcmdWithADouble("/adjoint_physics/SetCSBiasingFactorForCompton",this);
  SetCSBiasingFactorComptonCmd->SetGuidance("Set the CS biading factor for the adjoint Compton Cross section");
  SetCSBiasingFactorComptonCmd->SetParameterName("biasing_factor",false);
  SetCSBiasingFactorComptonCmd->AvailableForStates(G4State_PreInit);
  
  SetCSBiasingFactorBremCmd =  new G4UIcmdWithADouble("/adjoint_physics/SetCSBiasingFactorForBrem",this);
  SetCSBiasingFactorBremCmd->SetGuidance("Set the CS biading factor for the adjoint Brem Cross section");
  SetCSBiasingFactorBremCmd->SetParameterName("biasing_factor",false);
  SetCSBiasingFactorBremCmd->AvailableForStates(G4State_PreInit);
  
  SetCSBiasingFactorIonisationCmd =  new G4UIcmdWithADouble("/adjoint_physics/SetCSBiasingFactorForIonisation",this);
  SetCSBiasingFactorIonisationCmd->SetGuidance("Set the CS biading factor for the adjoint Ionisation Cross section");
  SetCSBiasingFactorIonisationCmd->SetParameterName("biasing_factor",false);
  SetCSBiasingFactorIonisationCmd->AvailableForStates(G4State_PreInit);
  
  SetCSBiasingFactorPEeffectCmd =  new G4UIcmdWithADouble("/adjoint_physics/SetCSBiasingFactorForPEeffect",this);
  SetCSBiasingFactorPEeffectCmd->SetGuidance("Set the CS biading factor for the adjoint PEeffect Cross section");
  SetCSBiasingFactorPEeffectCmd->SetParameterName("biasing_factor",false);
  SetCSBiasingFactorPEeffectCmd->AvailableForStates(G4State_PreInit);
   
 
#endif  
*/
  
 
  
  
    
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//

G4AdjointPhysicsMessenger::~G4AdjointPhysicsMessenger()
{

   
  
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//

void G4AdjointPhysicsMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  
  
  
  if ( command==UsepIonisationCmd){
        thePhysicsList->SetUseProtonIonisation(UsepIonisationCmd->GetNewBoolValue(newValue));
     
  }
  else if ( command==UseBremCmd){
  	thePhysicsList->SetUseBrem(UseBremCmd->GetNewBoolValue(newValue));
  }
  else if ( command==UseComptonCmd){
 	thePhysicsList->SetUseCompton(UseComptonCmd->GetNewBoolValue(newValue));
  }
  else if ( command==UseMSCmd){
    thePhysicsList->SetUseMS(UseMSCmd->GetNewBoolValue(newValue));
  }
  else if ( command==UsePEEffectCmd){
     thePhysicsList->SetUsePEEffect(UsePEEffectCmd->GetNewBoolValue(newValue));
  }
  else if ( command==UseGammaConversionCmd){
     thePhysicsList->SetUseGammaConversion(UseGammaConversionCmd->GetNewBoolValue(newValue));
  }
  else if ( command==UseEgainFluctuationCmd){
     thePhysicsList->SetUseEgainFluctuation(UseEgainFluctuationCmd->GetNewBoolValue(newValue));
  }

  else if ( command== SetEminAdjModelsCmd){
    thePhysicsList->SetEminAdjModels(SetEminAdjModelsCmd->GetNewDoubleValue(newValue));
  }
  else if ( command== SetEmaxAdjModelsCmd){
    thePhysicsList->SetEmaxAdjModels(SetEmaxAdjModelsCmd->GetNewDoubleValue(newValue));
  }
/*  
#ifdef TEST_MODE
  else if ( command== SetCSBiasingFactorComptonCmd){
    thePhysicsList->SetCSBiasingFactorCompton(SetCSBiasingFactorComptonCmd->GetNewDoubleValue(newValue));	
  }
  else if ( command== SetCSBiasingFactorBremCmd){
    thePhysicsList->SetCSBiasingFactorBrem(SetCSBiasingFactorBremCmd->GetNewDoubleValue(newValue));	
  }
  else if ( command== SetCSBiasingFactorIonisationCmd){
    thePhysicsList->SetCSBiasingFactorIonisation(SetCSBiasingFactorIonisationCmd->GetNewDoubleValue(newValue));	
  }
  else if ( command== SetCSBiasingFactorPEeffectCmd){
    thePhysicsList->SetCSBiasingFactorPEeffect(SetCSBiasingFactorPEeffectCmd->GetNewDoubleValue(newValue));	
  }
  
 
#endif    
*/  

}

