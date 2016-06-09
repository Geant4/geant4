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
// $Id: G4VisCommandsListManager.hh,v 1.1 2005/11/21 05:45:42 tinslay Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// Jane Tinslay, John Allison, Joseph Perl October 2005
//
// Class Description:
// Templated list manager commands which control list manager listing and
// selection.
// Class Description - End:

#ifndef G4VISCOMMANDLISTMANAGER_HH
#define G4VISCOMMANDLISTMANAGER_HH

#include "G4UImessenger.hh"
#include "G4String.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcommand.hh"

template <typename Manager>
class G4VisCommandListManagerList : public G4UImessenger {

public: // With description

  G4VisCommandListManagerList(Manager*, const G4String& placement);
  // Input list manager and command placement

  virtual ~G4VisCommandListManagerList();

  G4String GetCurrentValue(G4UIcommand*);
  void SetNewValue(G4UIcommand* command, G4String newValue);

  G4String Placement() const;

private:

  // Data members
  Manager* fpManager;
  G4String fPlacement;

  G4UIcmdWithAString* fpCommand;

};

// List command
template <typename Manager>
G4VisCommandListManagerList<Manager>::G4VisCommandListManagerList(Manager* manager, const G4String& placement)
  :fpManager(manager)
  ,fPlacement(placement)
{  
  G4String command = Placement()+"/list";
  
  fpCommand = new G4UIcmdWithAString(command, this);      
  fpCommand->SetGuidance("List objects registered with list manager");
  fpCommand->SetParameterName("name", true);       
}

template <typename Manager>
G4VisCommandListManagerList<Manager>::~G4VisCommandListManagerList()
{
  delete fpCommand;
}

template <typename Manager>
G4String
G4VisCommandListManagerList<Manager>::Placement() const
{
  return fPlacement;
}

template <typename Manager>
G4String 
G4VisCommandListManagerList<Manager>::GetCurrentValue(G4UIcommand*) 
{
  return "";
}

template <typename Manager>
void G4VisCommandListManagerList<Manager>::SetNewValue(G4UIcommand*, G4String name) 
{
  G4cout<<"Listing models available in "<<Placement()<<G4endl;

  assert (0 != fpManager);
  fpManager->Print(G4cout, name);
}    

//Select command
template <typename Manager>
class G4VisCommandListManagerSelect : public G4UImessenger {

public: // With description

  G4VisCommandListManagerSelect(Manager*, const G4String& placement);
  // Input list manager and command placement

  virtual ~G4VisCommandListManagerSelect();

  G4String GetCurrentValue(G4UIcommand*);
  void SetNewValue (G4UIcommand* command, G4String newValue);

private:

  Manager* fpManager;
  G4String fPlacement;

  G4UIcmdWithAString* fpCommand;

};

template <typename Manager>
G4VisCommandListManagerSelect<Manager>::G4VisCommandListManagerSelect(Manager* manager, const G4String& placement)
  :fpManager(manager)
  ,fPlacement(placement)
{  
  G4String command = placement+"/select";
  G4String guidance = "Select created object";
 
  fpCommand = new G4UIcmdWithAString(command, this);      
  fpCommand->SetGuidance(guidance);
  fpCommand->SetParameterName("name", false);       
}

template <typename Manager>
G4VisCommandListManagerSelect<Manager>::~G4VisCommandListManagerSelect()
{
  delete fpCommand;
}

template <typename Manager>
G4String 
G4VisCommandListManagerSelect<Manager>::GetCurrentValue(G4UIcommand*) 
{
  return "";
}

template <typename Manager>
void G4VisCommandListManagerSelect<Manager>::SetNewValue(G4UIcommand*, G4String name) 
{
  assert (0 != fpManager);
  fpManager->SetCurrent(name);
}    

#endif
