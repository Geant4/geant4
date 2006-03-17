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
// $Id: G4ModelCommandsT.hh,v 1.1 2006-03-17 03:24:02 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Generic model messenges providing set, setRGBA, setDefault and 
// setDefaultRGBA messengers. 
//
// Jane Tinslay March 2006
#ifndef G4MODELCOMMANDST_HH
#define G4MODELCOMMANDST_HH

#include "G4Colour.hh"
#include "G4String.hh"
#include "G4UIcommand.hh"
#include "G4VModelCommand.hh"

//Set variable colour
template <typename M>
class G4ModelCommandSet : public G4VModelCommand<M> {

public: // With description

  G4ModelCommandSet(M* model, const G4String& placement);

  virtual ~G4ModelCommandSet();

  void SetNewValue(G4UIcommand* command, G4String newValue);

private:

  G4UIcommand* stringCmd;
  G4UIcommand* componentCmd;

};

template <typename M>
G4ModelCommandSet<M>::G4ModelCommandSet(M* model, const G4String& placement)
  :G4VModelCommand<M>(model)
{
  //Set colour through a string
  G4String name = model->Name();
  G4String stringDir = placement+"/"+name+"/set";
  G4UIparameter* param(0);

  stringCmd = new G4UIcommand(stringDir, this);
  stringCmd->SetGuidance("Set colour through a string");   
  stringCmd->SetGuidance("Two inputs are expected.");
  
  param = new G4UIparameter("Variable ", 's', false);
  stringCmd->SetParameter(param);

  param = new G4UIparameter("Value ", 's', false);
  stringCmd->SetParameter(param);

  //Set colour through RGBA components
  G4String componentDir = placement+"/"+name+"/setRGBA";
  
  componentCmd = new G4UIcommand(componentDir, this);
  componentCmd->SetGuidance("Set colour through red, green, blue and alpha components");   
  componentCmd->SetGuidance("Five inputs are expected.");
  
  param = new G4UIparameter("Variable ", 's', false);
  componentCmd->SetParameter(param);

  param = new G4UIparameter("Red component ", 'd', false);
  componentCmd->SetParameter(param);

  param = new G4UIparameter("Green component ", 'd', false);
  componentCmd->SetParameter(param);

  param = new G4UIparameter("Blue component ", 'd', false);
  componentCmd->SetParameter(param);

  param = new G4UIparameter("Alpha component ", 'd', false);
  componentCmd->SetParameter(param);
} 

template <typename M>
G4ModelCommandSet<M>::~G4ModelCommandSet()
{ 
  delete stringCmd;
  delete componentCmd;
}

template <typename M>
void G4ModelCommandSet<M>::SetNewValue(G4UIcommand* cmd, G4String newValue)
{
  if (cmd == stringCmd) {
    G4String parameter;    
    G4String colour;
    std::istringstream is (newValue);
    is >> parameter >> colour;
    
    G4VModelCommand<M>::Model()->Set(parameter, colour);
  }

  if (cmd == componentCmd) {
    G4String parameter;    
    G4double red(0), green(0), blue(0), alpha(0);
    std::istringstream is (newValue);
    is >> parameter >> red >> green >> blue >> alpha;
    
    G4Colour myColour(red, green, blue, alpha);

    G4VModelCommand<M>::Model()->Set(parameter, myColour);
  }
}

//Set default colour
template <typename M>
class G4ModelCommandSetDefault : public G4VModelCommand<M> {

public: // With description

  G4ModelCommandSetDefault(M* model, const G4String& placement);

  virtual ~G4ModelCommandSetDefault();

  void SetNewValue(G4UIcommand* command, G4String newValue);

private:

  G4UIcommand* stringCmd;
  G4UIcommand* componentCmd;

};

template <typename M>
G4ModelCommandSetDefault<M>::G4ModelCommandSetDefault(M* model, const G4String& placement)
  :G4VModelCommand<M>(model)
{
  //Set colour through a string
  G4String name = model->Name();
  G4String stringDir = placement+"/"+name+"/setDefault";
  G4UIparameter* param(0);

  stringCmd = new G4UIcommand(stringDir, this);
  stringCmd->SetGuidance("Set default colour through a string");   
  stringCmd->SetGuidance("One input is expected.");

  param = new G4UIparameter("Value ", 's', false);
  stringCmd->SetParameter(param);

  //Set colour through RGBA components
  G4String componentDir = placement+"/"+name+"/setDefaultRGBA";
  
  componentCmd = new G4UIcommand(componentDir, this);
  componentCmd->SetGuidance("Set default colour through red, green, blue and alpha components.");   
  componentCmd->SetGuidance("Four inputs are expected.");

  param = new G4UIparameter("Red component ", 'd', false);
  componentCmd->SetParameter(param);

  param = new G4UIparameter("Green component ", 'd', false);
  componentCmd->SetParameter(param);

  param = new G4UIparameter("Blue component ", 'd', false);
  componentCmd->SetParameter(param);

  param = new G4UIparameter("Alpha component ", 'd', false);
  componentCmd->SetParameter(param);
} 

template <typename M>
G4ModelCommandSetDefault<M>::~G4ModelCommandSetDefault()
{ 
  delete stringCmd;
  delete componentCmd;
}

template <typename M>
void G4ModelCommandSetDefault<M>::SetNewValue(G4UIcommand* cmd, G4String newValue)
{
  if (cmd == stringCmd) {
    G4String colour;
    std::istringstream is (newValue);
    is >> colour;
    
    G4VModelCommand<M>::Model()->SetDefault(colour);
  }

  if (cmd == componentCmd) {
    G4double red(0), green(0), blue(0), alpha(0);
    std::istringstream is (newValue);
    is >> red >> green >> blue >> alpha;
    
    G4Colour myColour(red, green, blue, alpha);

    G4VModelCommand<M>::Model()->SetDefault(myColour);
  }
}

#endif
