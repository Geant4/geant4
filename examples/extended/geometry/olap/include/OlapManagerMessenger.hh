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
/// \file geometry/olap/include/OlapManagerMessenger.hh
/// \brief Definition of the OlapManagerMessenger class
//
//
// $Id$
//
// 
// --------------------------------------------------------------
// OlapManagerMessenger
//
// Author: Martin Liendl - Martin.Liendl@cern.ch
//
// --------------------------------------------------------------
//
#ifndef OlapManagerMessenger_h
#define OlapManagerMessenger_h

#include "G4UImessenger.hh"
#include "OlapLogManager.hh"

class OlapManager;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;
class G4UIcmdWith3Vector;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWithABool;

class OlapManagerMessenger : public G4UImessenger
{

public:
   OlapManagerMessenger(OlapManager*);
   ~OlapManagerMessenger();
   void SetNewValue(G4UIcommand*, G4String);
   
private:
   G4UIdirectory * theOlapDir;   
   G4UIcmdWith3VectorAndUnit * theRotationCmd;
   G4UIcmdWithoutParameter * theTriggerCmd;
   G4UIcmdWithAnInteger * theTriggerFullCmd;
   G4UIcmdWithAString * theGotoWorldCmd;
   G4UIcmdWithoutParameter * theLsCmd;
   G4UIcmdWithAString * theCdCmd;
   G4UIcmdWithAString * theListCmd;
   G4UIcmdWithoutParameter * thePwdCmd;
   G4UIcmdWithoutParameter * theWhereIsCmd;
   G4UIcmdWith3Vector * theSetGridCmd;
   G4UIcmdWithADoubleAndUnit * theDeltaCmd;
   G4UIcmdWithAString * theLogCmd;
   G4UIcmdWithAString * theLogByVolumeCmd;

   OlapManager * theManager;
   OlapLogManager * theLogManager;
   
};    
#endif
