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
// $Id: Tst23PhysicsListMessenger.cc,v 1.2 2002-12-05 02:19:06 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//---------------------------------------------------------------
//
//  Tst23PhysicsListMessenger.cc
//
//  Description:
//    This is a messenger class 
//
//---------------------------------------------------------------
#include "Tst23PhysicsListMessenger.hh"
#include "Tst23PhysicsList.hh"

#include "G4UImanager.hh"
#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"

#include "G4ios.hh"                 
#include "G4Tokenizer.hh"
#include "g4std/iomanip"               
#include "g4std/strstream"



//////////////////////////

Tst23PhysicsListMessenger::Tst23PhysicsListMessenger(Tst23PhysicsList* myPL)
         :myPhysList(myPL)
{
   myphysDir = new G4UIdirectory("/myphys/");
   myphysDir->SetGuidance("Physics setup commands.");

   cmdSetCut =  new G4UIcommand("/myphys/setCut",this);
   cmdSetCut->SetGuidance(" Set Cut In Range for region ");
   cmdSetCut->SetGuidance(" setCut region value unit");
   G4UIparameter * intParam = new G4UIparameter("region",'i',false);
   cmdSetCut->SetParameter(intParam);
   G4UIparameter * dblParam = new G4UIparameter("cut",'d',false);
   cmdSetCut->SetParameter(dblParam);
   G4UIparameter * untParam = new G4UIparameter("Unit",'s',false);
   cmdSetCut->SetParameter(untParam);
   cmdSetCut->AvailableForStates(G4State_PreInit,G4State_Idle);
}

Tst23PhysicsListMessenger::~Tst23PhysicsListMessenger()
{
  delete cmdSetCut;
  delete myphysDir;
}

void Tst23PhysicsListMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{
 if( command == cmdSetCut ){

  G4Tokenizer next( newValue );
   
   // check 1st argument
   G4int region;
   const char* t1 = G4String(next());
   G4std::istrstream is1((char*)t1);
   is1 >> region;
   
   // check 2nd argument
   G4double cut;
   const char* t2 = G4String(next());
   G4std::istrstream is2((char*)t2);
   is2 >> cut;

  // check 3rd argument
   G4String unit = G4String(next());

   cut *= G4UIcommand::ValueOf(unit);
 
   myPhysList->SetCutInRangeForRegion( cut, region);
 }  
}

G4String Tst23PhysicsListMessenger::GetCurrentValue(G4UIcommand * command)
{
  G4String cv;
  return cv;
}




