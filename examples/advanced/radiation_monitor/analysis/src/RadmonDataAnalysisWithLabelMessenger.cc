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
//[B
// File name:     RadmonDataAnalysisWithLabelMessenger.cc
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDataAnalysisWithLabelMessenger.cc,v 1.2 2006-06-28 13:44:41 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
//

// Messenger commands path
#define COMMANDS_PATH "/radmon/dataAnalysisFactory/"

// Include files
#include "RadmonDataAnalysisWithLabelMessenger.hh"



RadmonDataAnalysisWithLabelMessenger *        RadmonDataAnalysisWithLabelMessenger :: Instance(void)
{
 if (instance==0)
 {
  instance=new RadmonDataAnalysisWithLabelMessenger;
  
  if (!instance)
   G4Exception("RadmonDataAnalysisWithLabelMessenger::Instance: RadmonDataAnalysisWithLabelMessenger singleton not allocated.");
 }    
 
 return instance;
}




  
void                                            RadmonDataAnalysisWithLabelMessenger :: AddAvailableDataAnalysis(const G4String & name)
{
 availableDataAnalyses.insert(name);
}



void                                            RadmonDataAnalysisWithLabelMessenger :: RemoveAvailableDataAnalysis(const G4String & name)
{
 availableDataAnalyses.erase(name);
}





                                                RadmonDataAnalysisWithLabelMessenger :: RadmonDataAnalysisWithLabelMessenger()
:
 RadmonMessenger(COMMANDS_PATH, "Interactive labelled analysis data factory commands."),
 RADMON_INITIALIZE_COMMAND(Dump)
{
 RADMON_CREATE_COMMAND_0ARGS(Dump, "Print out the current analysis data constructors");
}



                                                RadmonDataAnalysisWithLabelMessenger :: ~RadmonDataAnalysisWithLabelMessenger()
{
 RADMON_DESTROY_COMMAND(Dump);
}





G4String                                        RadmonDataAnalysisWithLabelMessenger :: GetCurrentValue(G4UIcommand * /* command */)
{
 G4cout << "RadmonDataAnalysisWithLabelMessenger::GetCurrentValue(): Not supported" << G4endl;
 
 return G4String();
}



void                                            RadmonDataAnalysisWithLabelMessenger :: SetNewValue(G4UIcommand * command, G4String newValue)
{
 RADMON_BEGIN_LIST_SET_COMMANDS
  RADMON_SET_COMMAND(Dump)
 RADMON_END_LIST_SET_COMMANDS
}





// Events
void                                            RadmonDataAnalysisWithLabelMessenger :: OnDump(const G4String & /* value */)
{
 G4String indent("- ");
 
 G4String indent2(indent);
 indent2.prepend("  ");
 
 G4cout << indent << "Available analysis data:\n";
 
 AvailableDataAnalyses::const_iterator i(availableDataAnalyses.begin());
 AvailableDataAnalyses::const_iterator end(availableDataAnalyses.end());
 
 if (i==end)
  G4cout << indent2 << "No analysis data available.\n";
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





RadmonDataAnalysisWithLabelMessenger *        RadmonDataAnalysisWithLabelMessenger :: instance(0);
