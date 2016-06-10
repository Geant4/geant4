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
// $Id: G4VisCommandModelCreate.hh 66373 2012-12-18 09:41:34Z gcosmo $
//
// Jane Tinslay, John Allison, Joseph Perl October 2005
//
// Class Description:
// Templated create command for model factories. Factory generates new models
// and associated messengers. 
// Class Description - End:

#ifndef G4VISCOMMANDSMODELCREATE_HH
#define G4VISCOMMANDSMODELCREATE_HH

#include "G4VVisCommand.hh"
#include "G4String.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include <vector>

template <typename Factory>
class G4VisCommandModelCreate : public G4VVisCommand {

public: // With description

  G4VisCommandModelCreate(Factory*, const G4String& placement);
  // Input factory and command placement

  virtual ~G4VisCommandModelCreate();

  G4String GetCurrentValue(G4UIcommand*);
  void SetNewValue (G4UIcommand* command, G4String newValue);

  G4String Placement() const;

private:

  G4String NextName();

  // Data members
  Factory* fpFactory;
  G4String fPlacement;
  G4int fId;
  G4UIcmdWithAString* fpCommand;
  std::vector<G4UIcommand*>   fDirectoryList;

};

template <typename Factory>
G4VisCommandModelCreate<Factory>::G4VisCommandModelCreate(Factory* factory, const G4String& placement)
  :fpFactory(factory)
  ,fPlacement(placement)
  ,fId(0)
{  
  G4String factoryName = factory->Name();

  G4String command = Placement()+"/create/"+factoryName; 
  G4String guidance = "Create a "+factoryName+" model and associated messengers.";

  fpCommand = new G4UIcmdWithAString(command, this);      
  fpCommand->SetGuidance(guidance);
  fpCommand->SetGuidance("Generated model becomes current.");  
  fpCommand->SetParameterName("model-name", true);    
}

template <typename Factory>
G4VisCommandModelCreate<Factory>::~G4VisCommandModelCreate()
{
  delete fpCommand;
  
  unsigned i(0);
  for (i=0; i<fDirectoryList.size(); ++i) {
    delete fDirectoryList[i];
  }
}

template <typename Factory>
G4String
G4VisCommandModelCreate<Factory>::Placement() const
{
  return fPlacement;
}

template <typename Factory>
G4String
G4VisCommandModelCreate<Factory>::NextName()
{
  std::ostringstream oss;
  oss <<fpFactory->Name()<<"-" << fId++;
  return oss.str();
}

template <typename Factory>
G4String 
G4VisCommandModelCreate<Factory>::GetCurrentValue(G4UIcommand*) 
{
  return "";
}

template <typename Factory>
void G4VisCommandModelCreate<Factory>::SetNewValue(G4UIcommand*, G4String newName) 
{
  if (newName.isNull()) newName = NextName();

  assert (0 != fpFactory);

  // Create directory for new model commands
  G4String title = Placement()+"/"+newName+"/";
  G4String guidance = "Commands for "+newName+" model.";

  G4UIcommand* directory = new G4UIdirectory(title);
  directory->SetGuidance(guidance);
  fDirectoryList.push_back(directory);   

  // Create the model.
  typename Factory::ModelAndMessengers creation = fpFactory->Create(Placement(), newName);

  // Register model with vis manager
  fpVisManager->RegisterModel(creation.first);

  // Register associated messengers with vis manager
  typename Factory::Messengers::iterator iter = creation.second.begin();

  while (iter != creation.second.end()) {
    fpVisManager->RegisterMessenger(*iter);
    iter++;
  }
}    

#endif
