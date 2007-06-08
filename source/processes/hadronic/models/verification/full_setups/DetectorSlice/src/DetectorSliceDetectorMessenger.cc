#include "DetectorSliceDetectorMessenger.hh"

#include "DetectorSliceDetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "globals.hh"


DetectorSliceDetectorMessenger::
DetectorSliceDetectorMessenger(DetectorSliceDetectorConstruction* myDet)
  : theDetector(myDet) { 

  theDetectorDir = new G4UIdirectory("/mydet/");
  theDetectorDir->SetGuidance("Detector control.");
  
  theIsEmCalHomogeneous = new G4UIcmdWithABool("/mydet/isEmCalHomogeneous",this);
  theIsEmCalHomogeneous->SetParameterName("choiceIsEmCalHomogeneous",true);
  theIsEmCalHomogeneous->SetGuidance("Is the EM calorimeter homogeneous?");
  theIsEmCalHomogeneous->SetGuidance(" -> yes|y|true|t|1 : Homogeneous calorimeter ");
  theIsEmCalHomogeneous->SetGuidance(" -> no|n|false|f|0 : Sampling calorimeter ");
  theIsEmCalHomogeneous->SetDefaultValue( false ); // default: sampling calorimeter
  theIsEmCalHomogeneous->AvailableForStates(G4State_PreInit,G4State_Idle);  

  theIsHadCalHomogeneous = new G4UIcmdWithABool("/mydet/isHadCalHomogeneous",this);
  theIsHadCalHomogeneous->SetParameterName("choiceIsHadCalHomogeneous",true);
  theIsHadCalHomogeneous->SetGuidance("Is the HAD calorimeter homogeneous?");
  theIsHadCalHomogeneous->SetGuidance(" -> yes|y|true|t|1 : Homogeneous calorimeter ");
  theIsHadCalHomogeneous->SetGuidance(" -> no|n|false|f|0 : Sampling calorimeter ");
  theIsHadCalHomogeneous->SetDefaultValue( false ); // default: sampling calorimeter
  theIsHadCalHomogeneous->AvailableForStates(G4State_PreInit,G4State_Idle);  

  theTrackerMaterial = new G4UIcmdWithAString("/mydet/trackerMaterial",this);
  theTrackerMaterial->SetGuidance("Choice of the Tracker material:");
  theTrackerMaterial->SetGuidance("   silicon / scintillator / liquidArgon ");
  theTrackerMaterial->SetParameterName("choiceTrackerMaterial",true);
  theTrackerMaterial->SetDefaultValue("silicon");
  theTrackerMaterial->AvailableForStates(G4State_PreInit,G4State_Idle);

  theEmAbsorberMaterial = new G4UIcmdWithAString("/mydet/emAbsorberMaterial",this);
  theEmAbsorberMaterial->SetGuidance("Choice of the EM absorber material:");
  theEmAbsorberMaterial->SetGuidance("   iron / copper / tungsten / lead / PbWO4 / uranium ");
  theEmAbsorberMaterial->SetParameterName("choiceEmAbsorberMaterial",true);
  theEmAbsorberMaterial->SetDefaultValue("lead");
  theEmAbsorberMaterial->AvailableForStates(G4State_PreInit,G4State_Idle);

  theEmActiveMaterial = new G4UIcmdWithAString("/mydet/emActiveMaterial",this);
  theEmActiveMaterial->SetGuidance("Choice of the EM active material:");
  theEmActiveMaterial->SetGuidance("   scintillator / liquidArgon / PbWO4 / silicon / quartz ");
  theEmActiveMaterial->SetParameterName("choiceEmActiveMaterial",true);
  theEmActiveMaterial->SetDefaultValue("liquidArgon");
  theEmActiveMaterial->AvailableForStates(G4State_PreInit,G4State_Idle);

  theHadAbsorberMaterial = new G4UIcmdWithAString("/mydet/hadAbsorberMaterial",this);
  theHadAbsorberMaterial->SetGuidance("Choice of the HAD absorber material:");
  theHadAbsorberMaterial->SetGuidance("   iron / copper / tungsten / lead / PbWO4 / uranium ");
  theHadAbsorberMaterial->SetParameterName("choiceHadAbsorberMaterial",true);
  theHadAbsorberMaterial->SetDefaultValue("iron");
  theHadAbsorberMaterial->AvailableForStates(G4State_PreInit,G4State_Idle);

  theHadActiveMaterial = new G4UIcmdWithAString("/mydet/hadActiveMaterial",this);
  theHadActiveMaterial->SetGuidance("Choice of the HAD active material:");
  theHadActiveMaterial->SetGuidance("   scintillator / liquidArgon / PbWO4 / silicon / quartz ");
  theHadActiveMaterial->SetParameterName("choiceHadActiveMaterial",true);
  theHadActiveMaterial->SetDefaultValue("scintillator");
  theHadActiveMaterial->AvailableForStates(G4State_PreInit,G4State_Idle);

  theMuonMaterial = new G4UIcmdWithAString("/mydet/muonMaterial",this);
  theMuonMaterial->SetGuidance("Choice of the Muon material:");
  theMuonMaterial->SetGuidance("   iron / copper / tungsten / lead / uranium ");
  theMuonMaterial->SetParameterName("choiceMuonMaterial",true);
  theMuonMaterial->SetDefaultValue("iron");
  theMuonMaterial->AvailableForStates(G4State_PreInit,G4State_Idle);

  theTrackerLength = new G4UIcmdWithADouble("/mydet/trackerLength",this);
  theTrackerLength->SetParameterName("choiceTrackerLength",true);
  theTrackerLength->SetGuidance("Tracker length in [mm]");
  theTrackerLength->SetDefaultValue( 10.0 ); // default: 1 cm.
  theTrackerLength->AvailableForStates(G4State_PreInit,G4State_Idle);

  theEmAbsorberTotalLength = new G4UIcmdWithADouble("/mydet/emAbsorberTotalLength",this);
  theEmAbsorberTotalLength->SetParameterName("choiceEmAbsorberTotalLength",true);
  theEmAbsorberTotalLength->SetGuidance("EM Absorber total length, in [mm]");
  theEmAbsorberTotalLength->SetDefaultValue( 200.0 ); // default: 20 cm.
  theEmAbsorberTotalLength->AvailableForStates(G4State_PreInit,G4State_Idle);

  theHadAbsorberTotalLength = new G4UIcmdWithADouble("/mydet/hadAbsorberTotalLength",this);
  theHadAbsorberTotalLength->SetParameterName("choiceHadAbsorberTotalLength",true);
  theHadAbsorberTotalLength->SetGuidance("HAD Absorber total length, in [mm]");
  theHadAbsorberTotalLength->SetDefaultValue( 2000.0 ); // default: 2 meters.
  theHadAbsorberTotalLength->AvailableForStates(G4State_PreInit,G4State_Idle);

  theMuonLength = new G4UIcmdWithADouble("/mydet/muonLength",this);
  theMuonLength->SetParameterName("choiceMuonLength",true);
  theMuonLength->SetGuidance("Muon length, in [mm]");
  theMuonLength->SetDefaultValue( 1000.0 ); // default: 1 m.
  theMuonLength->AvailableForStates(G4State_PreInit,G4State_Idle);

  theDetectorRadius = new G4UIcmdWithADouble("/mydet/detectorRadius",this);
  theDetectorRadius->SetParameterName("choiceDetectorRadius",true);
  theDetectorRadius->SetGuidance("Detector Radius, in [mm]");
  theDetectorRadius->SetDefaultValue( 1000.0 ); // default: 1 m.
  theDetectorRadius->AvailableForStates(G4State_PreInit,G4State_Idle);

  theEmActiveLayerNumber = new G4UIcmdWithAnInteger("/mydet/emActiveLayerNumber",this);
  theEmActiveLayerNumber->SetParameterName("choiceEmActiveLayerNumber",true);
  theEmActiveLayerNumber->SetGuidance("Number of EM active layers");
  theEmActiveLayerNumber->SetDefaultValue( 50 );
  theEmActiveLayerNumber->AvailableForStates(G4State_PreInit,G4State_Idle);

  theHadActiveLayerNumber = new G4UIcmdWithAnInteger("/mydet/hadActiveLayerNumber",this);
  theHadActiveLayerNumber->SetParameterName("choiceHadActiveLayerNumber",true);
  theHadActiveLayerNumber->SetGuidance("Number of HAD active layers");
  theHadActiveLayerNumber->SetDefaultValue( 50 );
  theHadActiveLayerNumber->AvailableForStates(G4State_PreInit,G4State_Idle);

  theEmActiveLayerSize = new G4UIcmdWithADouble("/mydet/emActiveLayerSize",this);
  theEmActiveLayerSize->SetParameterName("choiceEmActiveLayerSize",true);
  theEmActiveLayerSize->SetGuidance("Size (thickness) of the EM active layer, in [mm]");
  theEmActiveLayerSize->SetDefaultValue( 1.0 ); // default: 1 mm.
  theEmActiveLayerSize->AvailableForStates(G4State_PreInit,G4State_Idle);

  theHadActiveLayerSize = new G4UIcmdWithADouble("/mydet/hadActiveLayerSize",this);
  theHadActiveLayerSize->SetParameterName("choiceHadActiveLayerSize",true);
  theHadActiveLayerSize->SetGuidance("Size (thickness) of the HAD active layer, in [mm]");
  theHadActiveLayerSize->SetDefaultValue( 4.0 ); // default: 4 millimeters.
  theHadActiveLayerSize->AvailableForStates(G4State_PreInit,G4State_Idle);

  theUpdateCommand = new G4UIcmdWithoutParameter("/mydet/update",this);
  theUpdateCommand->SetGuidance("Update calorimeter geometry.");
  theUpdateCommand->SetGuidance("This command MUST be applied before \"beamOn\" ");
  theUpdateCommand->SetGuidance("if you changed geometrical value(s).");
  theUpdateCommand->AvailableForStates(G4State_Idle);
 
}


DetectorSliceDetectorMessenger::~DetectorSliceDetectorMessenger() {

  delete theDetectorDir;

  delete theTrackerMaterial;
  delete theEmAbsorberMaterial;
  delete theEmActiveMaterial;
  delete theHadAbsorberMaterial;
  delete theHadActiveMaterial;
  delete theMuonMaterial;

  delete theIsEmCalHomogeneous;
  delete theIsHadCalHomogeneous;

  delete theTrackerLength;
  delete theEmAbsorberTotalLength;
  delete theHadAbsorberTotalLength;
  delete theMuonLength;
  delete theDetectorRadius;

  delete theEmActiveLayerNumber;
  delete theHadActiveLayerNumber;
  delete theEmActiveLayerSize;
  delete theHadActiveLayerSize;

  delete theUpdateCommand;

}


void DetectorSliceDetectorMessenger::
SetNewValue(G4UIcommand* command, G4String newValue) { 

  if ( command == theTrackerMaterial ) { 
    theDetector->SetTrackerMaterial( newValue );
  }

  if ( command == theEmAbsorberMaterial ) { 
    theDetector->SetEmAbsorberMaterial( newValue );
  }

  if ( command == theHadAbsorberMaterial ) { 
    theDetector->SetHadAbsorberMaterial( newValue );
  }

  if ( command == theEmActiveMaterial ) { 
    theDetector->SetEmActiveMaterial( newValue );
  }

  if ( command == theHadActiveMaterial ) { 
    theDetector->SetHadActiveMaterial( newValue );
  }

  if ( command == theMuonMaterial ) { 
    theDetector->SetMuonMaterial( newValue );
  }

  if ( command == theIsEmCalHomogeneous ) { 
    theDetector->SetIsEmCalHomogeneous( theIsEmCalHomogeneous->GetNewBoolValue(newValue) );
  }

  if ( command == theIsHadCalHomogeneous ) { 
    theDetector->SetIsHadCalHomogeneous( theIsHadCalHomogeneous->GetNewBoolValue(newValue) );
  }

  if ( command == theTrackerLength ) { 
    theDetector->
      SetTrackerLength( theTrackerLength->GetNewDoubleValue(newValue) );
  }

  if ( command == theEmAbsorberTotalLength ) { 
    theDetector->
      SetEmAbsorberTotalLength( theEmAbsorberTotalLength->GetNewDoubleValue(newValue) );
  }

  if ( command == theHadAbsorberTotalLength ) { 
    theDetector->
      SetHadAbsorberTotalLength( theHadAbsorberTotalLength->GetNewDoubleValue(newValue) );
  }

  if ( command == theMuonLength ) { 
    theDetector->
      SetMuonLength( theMuonLength->GetNewDoubleValue(newValue) );
  }

  if ( command == theDetectorRadius ) { 
    theDetector->
      SetDetectorRadius( theDetectorRadius->GetNewDoubleValue(newValue) );
  }

  if ( command == theEmActiveLayerNumber ) { 
    theDetector->SetEmActiveLayerNumber( theEmActiveLayerNumber->GetNewIntValue(newValue) );
  }

  if ( command == theHadActiveLayerNumber ) { 
    theDetector->SetHadActiveLayerNumber( theHadActiveLayerNumber->GetNewIntValue(newValue) );
  }

  if ( command == theEmActiveLayerSize ) { 
    theDetector->
      SetEmActiveLayerSize( theEmActiveLayerSize->GetNewDoubleValue(newValue) );
  }

  if ( command == theHadActiveLayerSize ) { 
    theDetector->
      SetHadActiveLayerSize( theHadActiveLayerSize->GetNewDoubleValue(newValue) );
  }

  if ( command == theUpdateCommand ) {
    theDetector->UpdateGeometry();
  }

}

