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
// $Id: OlapManagerMessenger.hh,v 1.1 2002-06-04 07:40:19 gcosmo Exp $
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
