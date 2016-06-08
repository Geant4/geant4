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
// $Id: B08MainMessenger.cc,v 1.1 2002/06/04 11:14:52 dressel Exp $
// GEANT4 tag $Name: geant4-04-01 $
//

#include "B08MainMessenger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcommand.hh"
#include "B08Scorer.hh"
#include "B08ScorePrinter.hh"
#include "G4VIStore.hh"

B08MainMessenger::
B08MainMessenger() : 
  fResultsFileName("test.txt"), 
  fNumberOfEvents(10),
  fUpper(1.),
  fLower(1.)
{
  
  fCmdFn = new G4UIcmdWithAString("/B08/resultsfile", this);
  fCmdFn->SetGuidance("set file name for output file");
  fCmdFn->SetParameterName("Messege",true);
  fCmdFn->SetDefaultValue("test.txt");

  fCmdEv = new G4UIcmdWithAnInteger("/B08/nevents", this);
  fCmdEv->SetGuidance("set number of events");
  fCmdEv->SetParameterName("Messege",true);
  fCmdEv->SetDefaultValue(10);

  fCmdUpper = new G4UIcmdWithADouble("/B08/upperlimit", this);
  fCmdUpper->SetGuidance("sets the upper limit for the weight window algorithm");
  fCmdUpper->SetParameterName("Messege",true);
  fCmdUpper->SetDefaultValue(1.);

  fCmdLower = new G4UIcmdWithADouble("/B08/lowerlimit", this);
  fCmdLower->SetGuidance("sets the lower limit for the weight window algorithm");
  fCmdLower->SetParameterName("Messege",true);
  fCmdLower->SetDefaultValue(1.);

}

void B08MainMessenger::SetNewValue(G4UIcommand* pCmd,G4String szValue){
  if(pCmd==fCmdFn){
    fResultsFileName = szValue;
  }
  if(pCmd==fCmdEv){
    fNumberOfEvents = atoi(szValue);
  }
  if(pCmd==fCmdUpper){
    fUpper = fCmdUpper->GetNewDoubleValue(szValue);
  }
  if(pCmd==fCmdLower){
    fLower = fCmdLower->GetNewDoubleValue(szValue);
  }
  
}


