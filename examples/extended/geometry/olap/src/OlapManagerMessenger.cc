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
//
// $Id: OlapManagerMessenger.cc,v 1.4 2010-08-24 07:57:14 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
// OlapManagerMessenger
//
// Author: Martin Liendl - Martin.Liendl@cern.ch
//
// --------------------------------------------------------------
//
#include "OlapManagerMessenger.hh"
#include "OlapManager.hh"
#include "OlapLogManager.hh"

#include "G4GeometryTolerance.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"

OlapManagerMessenger::OlapManagerMessenger(OlapManager* aManager)
  : theManager(aManager)
{

  theLogManager = OlapLogManager::GetOlapLogManager();

  G4double delta = G4GeometryTolerance::GetInstance()->GetAngularTolerance();

  theOlapDir = new G4UIdirectory("/olap/");
  theOlapDir->SetGuidance("Overlap detection facility");
  
  theRotationCmd = new G4UIcmdWith3VectorAndUnit("/olap/rotate",this);
  theRotationCmd->SetGuidance("rotate the new-world");
  theRotationCmd->AvailableForStates(G4State_Idle);
  theRotationCmd->SetUnitCategory("Angle");
  theRotationCmd->SetParameterName("rotAxisTheta", "rotAxisPhi", "rotAngle",true,true);

  theTriggerCmd = new G4UIcmdWithoutParameter("/olap/trigger",this);
  theTriggerCmd->SetGuidance("starts a single mother-daughters overlap detection.");
  theTriggerCmd->AvailableForStates(G4State_Idle);

  theTriggerFullCmd = new G4UIcmdWithAnInteger("/olap/triggerFull",this);
  theTriggerFullCmd->SetGuidance("starts a series of scans (only where mothers have");
  theTriggerFullCmd->SetGuidance("daughters. (-1 ... full detector)");
  theTriggerFullCmd->SetDefaultValue(1);
  theTriggerFullCmd->SetParameterName("nr",true);
  theTriggerFullCmd->SetRange("nr>=-1");
  theTriggerFullCmd->AvailableForStates(G4State_Idle);
  
  theDeltaCmd = new G4UIcmdWithADoubleAndUnit("/olap/delta",this);
  theDeltaCmd->SetGuidance("set boundary tolerance for overlaps in units of length");
  theDeltaCmd->SetDefaultValue(delta);
  theDeltaCmd->SetParameterName("delta",true);
  theDeltaCmd->SetRange("delta>=1.e-9"); // current G4-accuracy
  theDeltaCmd->SetUnitCategory("Length");
  theDeltaCmd->AvailableForStates(G4State_Idle);
  
  theSetGridCmd = new G4UIcmdWith3Vector("/olap/grid",this);
  theSetGridCmd->SetGuidance("set the grid for the generator (x-, y-, z- grid)");
  theSetGridCmd->SetDefaultValue(G4ThreeVector(3.,3.,3.));
  theSetGridCmd->SetParameterName("xGrid", "yGrid", "zGrid", true);
  theSetGridCmd->SetRange("xGrid>2. && yGrid >2. && zGrid >2.");
  theSetGridCmd->AvailableForStates(G4State_Idle);
    
  thePwdCmd = new G4UIcmdWithoutParameter("/olap/pwd",this);
  thePwdCmd->SetGuidance("show the position in the logical volume hierachy of the new world");
  thePwdCmd->AvailableForStates(G4State_Idle);
  
  theLsCmd = new G4UIcmdWithoutParameter("/olap/ls",this);
  theLsCmd->SetGuidance("lists the logical daughters of the current NewWorld");
  theLsCmd->AvailableForStates(G4State_Idle);
  
  theListCmd = new G4UIcmdWithAString("/olap/list",this);
  theListCmd->SetGuidance("lists all logical volumes which name matches regexp");
  theListCmd->AvailableForStates(G4State_Idle);
  
  theCdCmd = new G4UIcmdWithAString("/olap/cd",this);
  theCdCmd->SetGuidance("change to NewWorld like unix-cd");
  theCdCmd->AvailableForStates(G4State_Idle);
  
  theGotoWorldCmd = new G4UIcmdWithAString("/olap/goto",this);
  theGotoWorldCmd->SetGuidance("setting first logical vol matching regexp as NewWorld");
  theGotoWorldCmd->AvailableForStates(G4State_Idle);

  theLogCmd = new G4UIcmdWithAString("/olap/log", this);
  theLogCmd->SetGuidance("puts output into a single logfile");
  theLogCmd->SetDefaultValue("olap.log");
  theLogCmd->SetParameterName("name of logfile",true);
  theLogCmd->AvailableForStates(G4State_Idle);

  theLogByVolumeCmd = new G4UIcmdWithAString("/olap/logByVolume", this);
  theLogByVolumeCmd->SetGuidance("puts output into a logfile for each volume");
  theLogByVolumeCmd->SetDefaultValue("");
  theLogByVolumeCmd->SetParameterName("path of logfile(s)",true);
  theLogByVolumeCmd->AvailableForStates(G4State_Idle);

}


OlapManagerMessenger::~OlapManagerMessenger()
{
    delete theTriggerCmd;
    delete theTriggerFullCmd;
    delete theWhereIsCmd;
    delete theListCmd;
    delete theSetGridCmd;
    delete theGotoWorldCmd;
    delete theCdCmd;
    delete thePwdCmd;
    delete theLsCmd;
    delete theDeltaCmd;
    delete theOlapDir;
    delete theLogCmd;
    delete theLogByVolumeCmd;
    delete theRotationCmd;
}        


void OlapManagerMessenger::SetNewValue(G4UIcommand* aCmd, G4String aVal)
{
   
 if( aCmd == theTriggerCmd )
   { theManager->TriggerRun();}
 
 if( aCmd == theTriggerFullCmd )
   { theManager->TriggerFull(theTriggerFullCmd->GetNewIntValue(aVal)); }

 if( aCmd == theLsCmd )
   { theManager->LsLV(); }
   
 if( aCmd == thePwdCmd )
   { theManager->PwdLV(); }
   
 if( aCmd == theListCmd )
   { theManager->ListLV(aVal); }
 
 if( aCmd == theGotoWorldCmd )
   { theManager->GotoLV(aVal); }
   
 if( aCmd == theCdCmd )
   { theManager->ChangeLV(aVal); }  
   
 if ( aCmd == theSetGridCmd )
   { G4ThreeVector vec = theSetGridCmd->GetNew3VectorValue(aVal);
     theManager->SetGrid(G4int(vec[0]),G4int(vec[1]),G4int(vec[2]));
   } 
     
 if ( aCmd == theDeltaCmd )
   { 
     theManager->SetDelta(theDeltaCmd->GetNewDoubleValue(aVal)); 
     G4cout << "boundary tolerance set to " << theManager->Delta()/mm << "mm." << G4endl;
   }  
 if ( aCmd == theLogCmd )
   { theLogManager->Logging(aVal); }
   
 if (aCmd == theLogByVolumeCmd )
   { theLogManager->LogByVolume(aVal); }
 
 if (aCmd == theRotationCmd )
   { G4ThreeVector vec = theRotationCmd->GetNew3VectorValue(aVal);
     theManager->SetRotation(vec[0],vec[1],vec[2]); }  
}
