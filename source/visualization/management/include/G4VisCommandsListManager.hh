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
// $Id: G4VisCommandsListManager.hh 66373 2012-12-18 09:41:34Z gcosmo $
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

// Mode command
template <typename Manager>
class G4VisCommandManagerMode : public G4UImessenger {

public: // With description

  G4VisCommandManagerMode(Manager*, const G4String& placement);

  virtual ~G4VisCommandManagerMode();

  G4String GetCurrentValue(G4UIcommand*);
  void SetNewValue (G4UIcommand* command, G4String newValue);

private:

  Manager* fpManager;
  G4String fPlacement;

  G4UIcmdWithAString* fpCommand;

};
template <typename Manager>
G4VisCommandManagerMode<Manager>::G4VisCommandManagerMode(Manager* manager, const G4String& placement)
  :fpManager(manager)
  ,fPlacement(placement)
{  
  G4String command = fPlacement+"/mode";
  
  fpCommand = new G4UIcmdWithAString(command, this);      
  fpCommand->SetGuidance("Set mode of operation");
  fpCommand->SetParameterName("mode", false);       
  fpCommand->SetCandidates("soft hard");       
}

template <typename Manager>
G4VisCommandManagerMode<Manager>::~G4VisCommandManagerMode()
{
  delete fpCommand;
}

template <typename Manager>
G4String 
G4VisCommandManagerMode<Manager>::GetCurrentValue(G4UIcommand*) 
{
  return "";
}

template <typename Manager>
void G4VisCommandManagerMode<Manager>::SetNewValue(G4UIcommand*, G4String name) 
{
  assert (0 != fpManager);
  fpManager->SetMode(name);
}    

#endif
