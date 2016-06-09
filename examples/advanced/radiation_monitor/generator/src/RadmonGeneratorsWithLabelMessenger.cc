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
// File name:     RadmonGeneratorsWithLabelMessenger.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonGeneratorsWithLabelMessenger.cc,v 1.3 2006/06/29 16:16:39 gunter Exp $
// Tag:           $Name: geant4-08-02 $
//

// Messenger commands path
#define COMMANDS_PATH "/radmon/generatorsFactory/"

// Include files
#include "RadmonGeneratorsWithLabelMessenger.hh"



RadmonGeneratorsWithLabelMessenger *            RadmonGeneratorsWithLabelMessenger :: Instance(void)
{
 if (instance==0)
 {
  instance=new RadmonGeneratorsWithLabelMessenger;
  
  if (!instance)
   G4Exception("RadmonGeneratorsWithLabelMessenger::Instance: RadmonGeneratorsWithLabelMessenger singleton not allocated.");
 }    
 
 return instance;
}




  
void                                            RadmonGeneratorsWithLabelMessenger :: AddAvailableGenerator(const G4String & name)
{
 availableGenerators.insert(name);
}



void                                            RadmonGeneratorsWithLabelMessenger :: RemoveAvailableGenerator(const G4String & name)
{
 availableGenerators.erase(name);
}





                                                RadmonGeneratorsWithLabelMessenger :: RadmonGeneratorsWithLabelMessenger()
:
 RadmonMessenger(COMMANDS_PATH, "Interactive primary generators factory commands."),
 RADMON_INITIALIZE_COMMAND(Dump)
{
 RADMON_CREATE_COMMAND_0ARGS(Dump, "Print out the current available primary generators");
}



                                                RadmonGeneratorsWithLabelMessenger :: ~RadmonGeneratorsWithLabelMessenger()
{
 RADMON_DESTROY_COMMAND(Dump);
}





G4String                                        RadmonGeneratorsWithLabelMessenger :: GetCurrentValue(G4UIcommand * /* command */)
{
 G4cout << "RadmonGeneratorsWithLabelMessenger::GetCurrentValue(): Not supported" << G4endl;
 
 return G4String();
}



void                                            RadmonGeneratorsWithLabelMessenger :: SetNewValue(G4UIcommand * command, G4String newValue)
{
 RADMON_BEGIN_LIST_SET_COMMANDS
  RADMON_SET_COMMAND(Dump)
 RADMON_END_LIST_SET_COMMANDS
}





// Events
void                                            RadmonGeneratorsWithLabelMessenger :: OnDump(const G4String & /* value */)
{
 G4String indent("- ");
 
 G4String indent2(indent);
 indent2.prepend("  ");
 
 G4cout << indent << "Available generators:\n";
 
 AvailableGenerators::const_iterator i(availableGenerators.begin());
 AvailableGenerators::const_iterator end(availableGenerators.end());
 
 if (i==end)
  G4cout << indent2 << "No generators available.\n";
 else
 {
  while (i!=end)
  {
   G4cout << indent2 << (*i) << '\n';
   i++;
  }
 }
 
 G4cout << G4endl;
}





RadmonGeneratorsWithLabelMessenger *            RadmonGeneratorsWithLabelMessenger :: instance(0);
