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

#include "PassiveCarbonBeamLineMessenger.hh"
#include "PassiveCarbonBeamLine.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4SystemOfUnits.hh"

PassiveCarbonBeamLineMessenger::PassiveCarbonBeamLineMessenger(PassiveCarbonBeamLine* beamLine)
:passiveCarbon(beamLine)

{
    beamLineDir = new G4UIdirectory("/beamLine/");
    beamLineDir -> SetGuidance("set specification of ripple filter");

    RippleFilterDir = new G4UIdirectory("/beamLine/RippleFilter/");
    RippleFilterDir -> SetGuidance("set specification of ripple filter");

    finalCollimatorDir = new G4UIdirectory("/beamLine/FinalCollimator/");
    finalCollimatorDir -> SetGuidance("set specification of final collimator");

    RippleFilterMatCmd = new G4UIcmdWithAString("/beamLine/RippleFilter/RFMat",this);
    RippleFilterMatCmd -> SetGuidance("Set material of ripple filter");
    RippleFilterMatCmd -> SetParameterName("choice",false);
    RippleFilterMatCmd -> AvailableForStates(G4State_Idle);

    RippleFilterXPositionCmd = new G4UIcmdWithADoubleAndUnit("/beamLine/RippleFilter/position",this);
    RippleFilterXPositionCmd -> SetGuidance("Set position of ripple filter");
    RippleFilterXPositionCmd -> SetParameterName("Size",false);
    RippleFilterXPositionCmd -> SetDefaultUnit("mm");  
    RippleFilterXPositionCmd -> SetUnitCandidates("mm cm m");  
    RippleFilterXPositionCmd -> AvailableForStates(G4State_Idle);

    innerRadiusFinalCollimatorCmd = new G4UIcmdWithADoubleAndUnit("/beamLine/FinalCollimator/halfInnerRad",this);
    innerRadiusFinalCollimatorCmd -> SetGuidance("Set size of inner radius ( max 21.5 mm)");
    innerRadiusFinalCollimatorCmd -> SetParameterName("Size",false);
    innerRadiusFinalCollimatorCmd -> SetDefaultUnit("mm");  
    innerRadiusFinalCollimatorCmd -> SetUnitCandidates("mm cm m");  
    innerRadiusFinalCollimatorCmd -> AvailableForStates(G4State_Idle);

    }

PassiveCarbonBeamLineMessenger::~PassiveCarbonBeamLineMessenger()
{ 
   
    delete innerRadiusFinalCollimatorCmd;
    delete RippleFilterXPositionCmd;
    delete RippleFilterMatCmd;
    delete finalCollimatorDir;
    delete RippleFilterDir;  
    delete beamLineDir;
    
}


void PassiveCarbonBeamLineMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
    if( command == RippleFilterMatCmd )
    { passiveCarbon -> SetRippleFilterMaterial(newValue);
    }

    else if( command == RippleFilterXPositionCmd )
    { passiveCarbon -> SetRippleFilterXPosition
	(RippleFilterXPositionCmd -> GetNewDoubleValue(newValue));
    }

    else if( command == innerRadiusFinalCollimatorCmd )
    { passiveCarbon -> SetInnerRadiusFinalCollimator
	(innerRadiusFinalCollimatorCmd -> GetNewDoubleValue(newValue));
    }

    
}

