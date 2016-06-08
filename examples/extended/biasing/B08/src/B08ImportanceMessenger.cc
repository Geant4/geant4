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
// $Id: B08ImportanceMessenger.cc,v 1.1 2002/06/04 11:14:52 dressel Exp $
// GEANT4 tag $Name: geant4-04-01 $
//

#include "B08ImportanceMessenger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcommand.hh"
#include "B08Scorer.hh"
#include "B08ScorePrinter.hh"
#include "G4VIStore.hh"

#include "G4Pstring.hh"

B08ImportanceMessenger::
B08ImportanceMessenger(G4VIStore &aIStore) :
  fIStore(aIStore)
{}

void B08ImportanceMessenger::AddCell(const G4String &cellname, 
				     G4VPhysicalVolume *p){

  fCellVolMap[cellname] = p;
  
  G4String command = "/B08/importance/" + cellname + "/base";
  G4UIcmdWithADouble *base_cmd =
    new G4UIcmdWithADouble(command, this);

  
  fBaseMap[base_cmd] = p;
  
  command = "/B08/importance/" + cellname + "/exp";
  G4UIcmdWithADouble *expo_cmd =
    new G4UIcmdWithADouble(command, this);

  fExpMap[expo_cmd] = p;
  
  // create an element in the =base exponent map
  fBaseExpMap[p];
  
  fIStore.AddImportanceRegion(fBaseExpMap[p].GetImportance(), *p);

}


void B08ImportanceMessenger::
SetImportanceBase(const G4String &cellname, G4double base){
  G4VPhysicalVolume *p = fCellVolMap[cellname];
  B08BaseExp &be = fBaseExpMap[p];
  be.SetBase(base);
  fIStore.ChangeImportance(be.GetImportance(), *p);
}
void B08ImportanceMessenger::
SetImportanceExponent(const G4String &cellname, G4double exp){
  G4VPhysicalVolume *p = fCellVolMap[cellname];
  B08BaseExp &be = fBaseExpMap[p];
  be.SetExponent(exp);
  fIStore.ChangeImportance(be.GetImportance(), *p);
}


void B08ImportanceMessenger::SetNewValue(G4UIcommand* pCmd,G4String szValue){
  
  for (B08BaseMap::iterator itbase = fBaseMap.begin();
       itbase != fBaseMap.end(); itbase++){
    if (pCmd == itbase->first){
      B08BaseExpMap::iterator curent = fBaseExpMap.find(itbase->second);
      if (curent==fBaseExpMap.end()) {
	G4Exception("B08ImportanceMessenger: volume not found in BaseExpMap");
      }                            
      curent->second.SetBase(itbase->first->GetNewDoubleValue(szValue));
      fIStore.ChangeImportance(curent->second.GetImportance(), 
			       *(itbase->second));
      return;
    }
  }
  
  for (B08ExpMap::iterator itexp = fExpMap.begin();
       itexp != fExpMap.end(); itexp++){
    if (pCmd == itexp->first){
       B08BaseExpMap::iterator curent = fBaseExpMap.find(itexp->second);
       if (curent==fBaseExpMap.end()) {
	 G4Exception("B08ImportanceMessenger: volume not found in BaseExpMap");
       }
       curent->second.SetExponent(itexp->first->GetNewDoubleValue(szValue));
       fIStore.ChangeImportance(curent->second.GetImportance(), 
				*(itexp->second));
       return;
    }
  }
  
}


