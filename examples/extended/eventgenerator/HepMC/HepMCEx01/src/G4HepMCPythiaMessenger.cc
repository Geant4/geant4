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

// ====================================================================
//
//   G4HepMCPythiaMessenger.cc
//   $Id: G4HepMCPythiaMessenger.cc,v 1.1 2002-04-29 20:44:15 asaim Exp $
//
// ====================================================================
#include "G4HepMCPythiaMessenger.hh"
#include "G4HepMCPythiaInterface.hh"

#include "g4std/fstream"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"

////////////////////////////////////////////////////////////////////////////
G4HepMCPythiaMessenger::G4HepMCPythiaMessenger(G4HepMCPythiaInterface* agen)
  : gen(agen)
////////////////////////////////////////////////////////////////////////////
{
  dir= new G4UIdirectory("/generator/pythia/");
  dir-> SetGuidance("Commands for Pythia event generation");

  verbose= new G4UIcmdWithAnInteger("/generator/pythia/verbose",this);
  verbose-> SetGuidance("set verbose level");
  verbose-> SetParameterName("verboseLevel", false, false);
  verbose-> SetRange("verboseLevel>=0 && verboseLevel<=2");

  mpylist= new G4UIcmdWithAnInteger("/generator/pythia/pylist",this);
  mpylist-> SetGuidance("set argument of pylist (not called if mlist=0)");
  mpylist-> SetParameterName("mlist", false, false);
  mpylist-> SetRange("mlist>=0 && mlist<=3");

  print= new G4UIcmdWithoutParameter("/generator/pythia/print", this);
  print-> SetGuidance("print user information.");

  cpyinit= new G4UIcommand("/generator/pythia/pyinit", this);
  cpyinit-> SetGuidance("call PYINIT");
  G4UIparameter* frame=
    new G4UIparameter("frame of the experiment", 's', false);
  cpyinit-> SetParameter(frame);
  G4UIparameter* beam= new G4UIparameter("beam particle", 's', false);
  cpyinit-> SetParameter(beam);
  G4UIparameter* target= new G4UIparameter("target particle", 's', false);
  cpyinit-> SetParameter(target);
  G4UIparameter* win= new G4UIparameter("energy of system (GeV)", 'd', false);
  cpyinit-> SetParameter(win);

  cpystat= new G4UIcmdWithAnInteger("/generator/pythia/pystat", this);
  cpystat-> SetGuidance("call PYSTAT");
  cpystat-> SetParameterName("mstat", false, false);
  cpystat-> SetRange("mstat>=1 && mstat<=5");

  cpygive= new G4UIcommand("/generator/pythia/pygive",this);
  cpygive-> SetGuidance("call PYGIVE");
  G4UIparameter* parameter= new G4UIparameter ("Parameter", 's', false);
  cpygive-> SetParameter(parameter);  

  setUserParameters= 
    new G4UIcmdWithoutParameter("/generator/pythia/setUserParameters",this);
  setUserParameters-> 
    SetGuidance("Set user parameters in the Pythia common blocks");

  setSeed= new G4UIcmdWithAnInteger("/generator/pythia/setSeed", this);
  setSeed-> SetGuidance("set initial seed.");

  cpyrget= new G4UIcommand("/generator/pythia/pyrget", this);
  cpyrget-> SetGuidance("call PYRGET");
  G4UIparameter* lun, *move;
  lun= new G4UIparameter("logical file number", 'i', false);
  cpyrget-> SetParameter(lun);
  move= new G4UIparameter("choice of adding a new record", 'i', true);
  move-> SetDefaultValue(-1);
  cpyrget-> SetParameter(move);

  cpyrset= new G4UIcommand("/generator/pythia/pyrset", this);
  cpyrset-> SetGuidance("call PYRSET");
  lun= new G4UIparameter("logical file number", 'i', false);
  cpyrset-> SetParameter(lun);
  move= new G4UIparameter("choice of adding a new record", 'i', true);
  move-> SetDefaultValue(0);
  cpyrset-> SetParameter(move);

  printRandomStatus= 
    new G4UIcmdWithAString("/generator/pythia/printRandomStatus", this);
  printRandomStatus-> SetGuidance("print random number status.");
  printRandomStatus-> SetParameterName("filename", true, false);
  printRandomStatus-> SetDefaultValue("std::cout");
}

/////////////////////////////////////////////////
G4HepMCPythiaMessenger::~G4HepMCPythiaMessenger()
/////////////////////////////////////////////////
{
  delete verbose;
  delete mpylist;
  delete print;
  delete cpyinit;
  delete cpystat;  
  delete cpygive;
  delete setUserParameters;
  delete setSeed;
  delete cpyrget;
  delete cpyrset;
  delete printRandomStatus;

  delete dir;
}

//////////////////////////////////////////////////////////////
void G4HepMCPythiaMessenger::SetNewValue(G4UIcommand* command, 
					 G4String newValues)
//////////////////////////////////////////////////////////////
{
  if(command == verbose) {  // /verbose ...
    G4int level= verbose-> GetNewIntValue(newValues);
    gen-> SetVerboseLevel(level);

  } else if (command == mpylist) { // /mpylist ...
    G4int mlist= mpylist-> GetNewIntValue(newValues);
    gen-> SetPylist(mlist);

  } else if (command == print) { // /print ...
    gen-> Print();

  } else if (command == cpyinit) { // /pyinit ...
    const char* strvaluelist= newValues.c_str();
    G4std::istrstream is((char*)strvaluelist);
    G4String sframe, sbeam, starget; G4double dwin;
    is >> sframe >> sbeam >> starget >> dwin;
    gen-> CallPyinit(sframe, sbeam, starget, dwin);

  } else if (command == cpystat) { // /pystat ...
    G4int imod= cpystat-> GetNewIntValue(newValues);
    gen-> CallPystat(imod);

  } else if (command == cpygive) { // /pygive ...
    G4String s= newValues;
    gen-> CallPygive(s);

  } else if (command == setUserParameters) { // /setUserParameters ...
    gen-> SetUserParameters();

  } else if (command == setSeed) { // /setSeed ...
    G4int iseed= setSeed-> GetNewIntValue(newValues);
    gen-> SetRandomSeed(iseed);

  } else if (command == cpyrget) { // /pyrget ...
    const char* strvaluelist= newValues.c_str();
    G4std::istrstream is((char*)strvaluelist);
    G4int lun, move;
    is >> lun >> move;
    gen-> CallPyrget(lun, move);

  } else if (command == cpyrset) { // /pyrset ...
    const char* strvaluelist= newValues.c_str();
    G4std::istrstream is((char*)strvaluelist);
    G4int lun, move;
    is >> lun >> move;
    gen-> CallPyrset(lun, move);

  } else if (command == printRandomStatus) { // /printRandomStatus ...
    G4String s= newValues;
    if (newValues == "std::cout") {
      gen-> PrintRandomStatus();
    } else {
      // to a file (overwrite mode)
      G4std::ofstream ofs;
      ofs.open(s.c_str(), std::ios::out);  
      //ofs.open(randomStatusFileName.c_str(), std::ios::out|std::ios::app);
      ofs.setf(std::ios::fixed | std::ios::showpoint);
      gen-> PrintRandomStatus(ofs);
      ofs.close();
    }
  }
}

//////////////////////////////////////////////////////////////////////
G4String G4HepMCPythiaMessenger::GetCurrentValue(G4UIcommand* command)
//////////////////////////////////////////////////////////////////////
{
  G4String cv;
  if (command == verbose) {
    cv= verbose-> ConvertToString(gen->GetVerboseLevel());
  } else  if (command == mpylist) {
    cv= verbose-> ConvertToString(gen->GetPylist());
  }
  return cv;
}

