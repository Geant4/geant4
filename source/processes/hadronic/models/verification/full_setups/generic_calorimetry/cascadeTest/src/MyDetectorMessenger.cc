#include "MyDetectorMessenger.hh"

#include "MyDetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAString.hh"
#include "globals.hh"


MyDetectorMessenger::MyDetectorMessenger(MyDetectorConstruction* myDet)
  : theDetector(myDet) { 

  theDetectorDir = new G4UIdirectory("/mydet/");
  theDetectorDir->SetGuidance("Detector control.");
  
  theFieldCommand = new G4UIcmdWithADoubleAndUnit("/mydet/setField",this);  
  theFieldCommand->SetGuidance("Define magnetic field.");
  theFieldCommand->SetGuidance("Magnetic field will be in X direction.");
  theFieldCommand->SetParameterName("Bx",false);
  theFieldCommand->SetDefaultUnit("tesla");
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

}


MyDetectorMessenger::~MyDetectorMessenger() {

  delete theFieldCommand;
  delete theDetectorDir;
  delete theAbsorberMaterial;
  delete theActiveMaterial;

}


void MyDetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue) { 

  if ( command == theFieldCommand ) { 
    theDetector->SetMagField( theFieldCommand->GetNewDoubleValue(newValue) );
  }

  if ( command == theAbsorberMaterial ) { 
    theDetector->SetAbsorberMaterial( newValue );
  }

  if ( command == theActiveMaterial ) { 
    theDetector->SetActiveMaterial( newValue );
  }

}

