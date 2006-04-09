#include "StatAccepTestDetectorMessenger.hh"

#include "StatAccepTestDetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "globals.hh"


StatAccepTestDetectorMessenger::
StatAccepTestDetectorMessenger(StatAccepTestDetectorConstruction* myDet)
  : theDetector(myDet) { 

  theDetectorDir = new G4UIdirectory("/mydet/");
  theDetectorDir->SetGuidance("Detector control.");
  
  theFieldCommand = new G4UIcmdWithADoubleAndUnit("/mydet/setField",this);  
  theFieldCommand->SetGuidance("Define uniform magnetic field along Y.");
  theFieldCommand->SetGuidance(" -> in unit of  [Tesla]");
  theFieldCommand->SetParameterName("By",false);
  theFieldCommand->SetDefaultValue( 0.0 );
  theFieldCommand->SetUnitCategory("Magnetic flux density");
  theFieldCommand->AvailableForStates(G4State_PreInit,G4State_Idle);  

  theAbsorberMaterial = new G4UIcmdWithAString("/mydet/absorberMaterial",this);
  theAbsorberMaterial->SetGuidance("Choice of the absorber material:");
  theAbsorberMaterial->SetGuidance("   iron / copper / tungsten / lead / PbWO4 / uranium ");
  theAbsorberMaterial->SetParameterName("choiceAbsorberMaterial",true);
  theAbsorberMaterial->SetDefaultValue("iron");
  theAbsorberMaterial->AvailableForStates(G4State_PreInit,G4State_Idle);

  theActiveMaterial = new G4UIcmdWithAString("/mydet/activeMaterial",this);
  theActiveMaterial->SetGuidance("Choice of the active material:");
  theActiveMaterial->SetGuidance("   scintillator / liquidArgon / PbWO4 / silicon / quartz ");
  theActiveMaterial->SetParameterName("choiceActiveMaterial",true);
  theActiveMaterial->SetDefaultValue("scintillator");
  theActiveMaterial->AvailableForStates(G4State_PreInit,G4State_Idle);

  theIsCalHomogeneous = new G4UIcmdWithABool("/mydet/isCalHomogeneous",this);
  theIsCalHomogeneous->SetParameterName("choiceIsCalHomogeneous",true);
  theIsCalHomogeneous->SetGuidance("Is the calorimeter homogeneous?");
  theIsCalHomogeneous->SetGuidance(" -> yes|y|true|t|1 : Homogeneous calorimeter ");
  theIsCalHomogeneous->SetGuidance(" -> no|n|false|f|0 : Sampling calorimeter ");
  theIsCalHomogeneous->SetDefaultValue( false ); // default: sampling calorimeter
  theIsCalHomogeneous->AvailableForStates(G4State_PreInit,G4State_Idle);  
  
  theIsUnitInLambda = new G4UIcmdWithABool("/mydet/isUnitInLambda",this);
  theIsUnitInLambda->SetParameterName("choiceIsUnitInLambda",true);
  theIsUnitInLambda->SetGuidance("Is unit for absorber length in lambda?");
  theIsUnitInLambda->SetGuidance(" -> yes|y|true|t|1 : unit in lambda");
  theIsUnitInLambda->SetGuidance(" -> no|n|false|f|0 : unit in [mm]");
  theIsUnitInLambda->SetDefaultValue( false ); // default: unit in [mm].
  theIsUnitInLambda->AvailableForStates(G4State_PreInit,G4State_Idle);  
  
  theAbsorberTotalLength = new G4UIcmdWithADouble("/mydet/absorberTotalLength",this);
  theAbsorberTotalLength->SetParameterName("choiceAbsorberTotalLength",true);
  theAbsorberTotalLength->SetGuidance("Absorber total length");
  theAbsorberTotalLength->SetGuidance(" -> in unit of lambda or [mm]");
  theAbsorberTotalLength->SetGuidance(" -> depending on value of choiceIsUnitInLambda");
  theAbsorberTotalLength->SetDefaultValue( 2000.0 ); // default: 2 meters.
  theAbsorberTotalLength->AvailableForStates(G4State_PreInit,G4State_Idle);

  theCalorimeterRadius = new G4UIcmdWithADouble("/mydet/calorimeterRadius",this);
  theCalorimeterRadius->SetParameterName("choiceCalorimeterRadius",true);
  theCalorimeterRadius->SetGuidance("Calorimeter Radius");
  theCalorimeterRadius->SetGuidance(" -> in unit of lambda or [mm]");
  theCalorimeterRadius->SetGuidance(" -> depending on value of choiceIsUnitInLambda");
  theCalorimeterRadius->SetDefaultValue( 1000.0 ); // default: 1 meter.
  theCalorimeterRadius->AvailableForStates(G4State_PreInit,G4State_Idle);

  theActiveLayerNumber = new G4UIcmdWithAnInteger("/mydet/activeLayerNumber",this);
  theActiveLayerNumber->SetParameterName("choiceActiveLayerNumber",true);
  theActiveLayerNumber->SetGuidance("Number of active layers");
  theActiveLayerNumber->SetDefaultValue( 50 );
  theActiveLayerNumber->AvailableForStates(G4State_PreInit,G4State_Idle);

  theActiveLayerSize = new G4UIcmdWithADouble("/mydet/activeLayerSize",this);
  theActiveLayerSize->SetParameterName("choiceActiveLayerSize",true);
  theActiveLayerSize->SetGuidance("Size (thickness) of the active layer, in [mm]");
  theActiveLayerSize->SetDefaultValue( 4.0 ); // default: 4 millimeters.
  theActiveLayerSize->AvailableForStates(G4State_PreInit,G4State_Idle);

  theReadoutLayerNumber = new G4UIcmdWithAnInteger("/mydet/readoutLayerNumber",this);
  theReadoutLayerNumber->SetParameterName("choiceReadoutLayerNumber",true);
  theReadoutLayerNumber->SetGuidance("Number of readout layers: it must be a");
  theReadoutLayerNumber->SetGuidance("divisor of the number of active layers.");
  theReadoutLayerNumber->SetDefaultValue( 50 );
  theReadoutLayerNumber->AvailableForStates(G4State_PreInit,G4State_Idle);

  theIsRadiusUnitInLambda = new G4UIcmdWithABool("/mydet/isRadiusUnitInLambda",this);
  theIsRadiusUnitInLambda->SetParameterName("choiceIsRadiusUnitInLambda",true);
  theIsRadiusUnitInLambda->SetGuidance("Is unit of radius in lambda?");
  theIsRadiusUnitInLambda->SetGuidance(" -> yes|y|true|t|1 : unit in lambda");
  theIsRadiusUnitInLambda->SetGuidance(" -> no|n|false|f|0 : unit in [mm]");
  theIsRadiusUnitInLambda->SetDefaultValue( false ); // default: unit in [mm].
  theIsRadiusUnitInLambda->AvailableForStates(G4State_PreInit,G4State_Idle);  
  
  theRadiusBinSize = new G4UIcmdWithADouble("/mydet/radiusBinSize",this);
  theRadiusBinSize->SetParameterName("choiceRadiusBinSize",true);
  theRadiusBinSize->SetGuidance("Size of the radius bin");
  theRadiusBinSize->SetGuidance(" -> in unit of lambda or [mm]");
  theRadiusBinSize->SetGuidance(" -> depending on value of choiceIsRadiusUnitInLambda");
  theRadiusBinSize->SetDefaultValue( 100.0 ); // default: 10 centimeters.
  theRadiusBinSize->AvailableForStates(G4State_PreInit,G4State_Idle);

  theRadiusBinNumber = new G4UIcmdWithAnInteger("/mydet/radiusBinNumber",this);
  theRadiusBinNumber->SetParameterName("choiceRadiusBinNumber",true);
  theRadiusBinNumber->SetGuidance("Number of radius bins");
  theRadiusBinNumber->SetGuidance(" NB) add one more to collect all hits with radius");
  theRadiusBinNumber->SetGuidance("     above the last (normal) one.");
  theRadiusBinNumber->SetDefaultValue( 11 );
  theRadiusBinNumber->AvailableForStates(G4State_PreInit,G4State_Idle);

  theUpdateCommand = new G4UIcmdWithoutParameter("/mydet/update",this);
  theUpdateCommand->SetGuidance("Update calorimeter geometry.");
  theUpdateCommand->SetGuidance("This command MUST be applied before \"beamOn\" ");
  theUpdateCommand->SetGuidance("if you changed geometrical value(s).");
  theUpdateCommand->AvailableForStates(G4State_Idle);
 
}


StatAccepTestDetectorMessenger::~StatAccepTestDetectorMessenger() {

  delete theFieldCommand;
  delete theDetectorDir;
  delete theAbsorberMaterial;
  delete theActiveMaterial;

  delete theIsCalHomogeneous;
  delete theIsUnitInLambda;
  delete theAbsorberTotalLength;
  delete theCalorimeterRadius;
  delete theActiveLayerNumber;
  delete theActiveLayerSize;
  delete theReadoutLayerNumber;

  delete theIsRadiusUnitInLambda;
  delete theRadiusBinSize;
  delete theRadiusBinNumber;

  delete theUpdateCommand;

}


void StatAccepTestDetectorMessenger::
SetNewValue(G4UIcommand* command, G4String newValue) { 

  if ( command == theFieldCommand ) { 
    theDetector->SetMagField( theFieldCommand->GetNewDoubleValue(newValue) );
  }

  if ( command == theAbsorberMaterial ) { 
    theDetector->SetAbsorberMaterial( newValue );
  }

  if ( command == theActiveMaterial ) { 
    theDetector->SetActiveMaterial( newValue );
  }

  if ( command == theIsCalHomogeneous ) { 
    theDetector->SetIsCalHomogeneous( theIsCalHomogeneous->GetNewBoolValue(newValue) );
  }

  if ( command == theIsUnitInLambda ) { 
    theDetector->SetIsUnitInLambda( theIsUnitInLambda->GetNewBoolValue(newValue) );
  }

  if ( command == theAbsorberTotalLength ) { 
    theDetector->
      SetAbsorberTotalLength( theAbsorberTotalLength->GetNewDoubleValue(newValue) );
  }

  if ( command == theCalorimeterRadius ) { 
    theDetector->
      SetCalorimeterRadius( theCalorimeterRadius->GetNewDoubleValue(newValue) );
  }

  if ( command == theActiveLayerNumber ) { 
    theDetector->SetActiveLayerNumber( theActiveLayerNumber->GetNewIntValue(newValue) );
  }

  if ( command == theActiveLayerSize ) { 
    theDetector->
      SetActiveLayerSize( theActiveLayerSize->GetNewDoubleValue(newValue) );
  }

  if ( command == theReadoutLayerNumber ) { 
    theDetector->SetReadoutLayerNumber( theReadoutLayerNumber->GetNewIntValue(newValue) );
  }

  if ( command == theIsRadiusUnitInLambda ) { 
    theDetector->
      SetIsRadiusUnitInLambda( theIsRadiusUnitInLambda->GetNewBoolValue(newValue) );
  }

  if ( command == theRadiusBinSize ) { 
    theDetector->SetRadiusBinSize( theRadiusBinSize->GetNewDoubleValue(newValue) );
  }
  
  if ( command == theRadiusBinNumber ) { 
    theDetector->SetRadiusBinNumber( theRadiusBinNumber->GetNewIntValue(newValue) );
  }

  if ( command == theUpdateCommand ) {
    theDetector->UpdateGeometry();
  }

}

