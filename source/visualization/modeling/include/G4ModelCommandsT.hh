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
// $Id: G4ModelCommandsT.hh,v 1.2 2006-03-28 18:01:18 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Generic model messenges. 
//
// Jane Tinslay March 2006
//
#ifndef G4MODELCOMMANDST_HH
#define G4MODELCOMMANDST_HH

#include "G4Colour.hh"
#include "G4String.hh"
#include "G4UIcmdWithABool.hh"
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

  G4UIcommand* fpStringCmd;
  G4UIcommand* fpComponentCmd;

};

template <typename M>
G4ModelCommandSet<M>::G4ModelCommandSet(M* model, const G4String& placement)
  :G4VModelCommand<M>(model)
{
  //Set colour through a string
  G4String name = model->Name();
  G4String stringDir = placement+"/"+name+"/set";
  G4UIparameter* param(0);

  fpStringCmd = new G4UIcommand(stringDir, this);
  fpStringCmd->SetGuidance("Set colour through a string");   
  fpStringCmd->SetGuidance("Two inputs are expected.");
  
  param = new G4UIparameter("Variable ", 's', false);
  fpStringCmd->SetParameter(param);

  param = new G4UIparameter("Value ", 's', false);
  fpStringCmd->SetParameter(param);

  //Set colour through RGBA components
  G4String componentDir = placement+"/"+name+"/setRGBA";
  
  fpComponentCmd = new G4UIcommand(componentDir, this);
  fpComponentCmd->SetGuidance("Set colour through red, green, blue and alpha components");   
  fpComponentCmd->SetGuidance("Five inputs are expected.");
  
  param = new G4UIparameter("Variable ", 's', false);
  fpComponentCmd->SetParameter(param);

  param = new G4UIparameter("Red component ", 'd', false);
  fpComponentCmd->SetParameter(param);

  param = new G4UIparameter("Green component ", 'd', false);
  fpComponentCmd->SetParameter(param);

  param = new G4UIparameter("Blue component ", 'd', false);
  fpComponentCmd->SetParameter(param);

  param = new G4UIparameter("Alpha component ", 'd', false);
  fpComponentCmd->SetParameter(param);
} 

template <typename M>
G4ModelCommandSet<M>::~G4ModelCommandSet()
{ 
  delete fpStringCmd;
  delete fpComponentCmd;
}

template <typename M>
void G4ModelCommandSet<M>::SetNewValue(G4UIcommand* cmd, G4String newValue)
{
  if (cmd == fpStringCmd) {
    G4String parameter;    
    G4String colour;
    std::istringstream is (newValue);
    is >> parameter >> colour;
    
    G4VModelCommand<M>::Model()->Set(parameter, colour);
  }

  if (cmd == fpComponentCmd) {
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

  G4UIcommand* fpStringCmd;
  G4UIcommand* fpComponentCmd;

};

template <typename M>
G4ModelCommandSetDefault<M>::G4ModelCommandSetDefault(M* model, const G4String& placement)
  :G4VModelCommand<M>(model)
{
  //Set colour through a string
  G4String name = model->Name();
  G4String stringDir = placement+"/"+name+"/setDefault";
  G4UIparameter* param(0);

  fpStringCmd = new G4UIcommand(stringDir, this);
  fpStringCmd->SetGuidance("Set default colour through a string");   
  fpStringCmd->SetGuidance("One input is expected.");

  param = new G4UIparameter("Value ", 's', false);
  fpStringCmd->SetParameter(param);

  //Set colour through RGBA components
  G4String componentDir = placement+"/"+name+"/setDefaultRGBA";
  
  fpComponentCmd = new G4UIcommand(componentDir, this);
  fpComponentCmd->SetGuidance("Set default colour through red, green, blue and alpha components.");   
  fpComponentCmd->SetGuidance("Four inputs are expected.");

  param = new G4UIparameter("Red component ", 'd', false);
  fpComponentCmd->SetParameter(param);

  param = new G4UIparameter("Green component ", 'd', false);
  fpComponentCmd->SetParameter(param);

  param = new G4UIparameter("Blue component ", 'd', false);
  fpComponentCmd->SetParameter(param);

  param = new G4UIparameter("Alpha component ", 'd', false);
  fpComponentCmd->SetParameter(param);
} 

template <typename M>
G4ModelCommandSetDefault<M>::~G4ModelCommandSetDefault()
{ 
  delete fpStringCmd;
  delete fpComponentCmd;
}

template <typename M>
void G4ModelCommandSetDefault<M>::SetNewValue(G4UIcommand* cmd, G4String newValue)
{
  if (cmd == fpStringCmd) {
    G4String colour;
    std::istringstream is (newValue);
    is >> colour;
    
    G4VModelCommand<M>::Model()->SetDefault(colour);
  }

  if (cmd == fpComponentCmd) {
    G4double red(0), green(0), blue(0), alpha(0);
    std::istringstream is (newValue);
    is >> red >> green >> blue >> alpha;
    
    G4Colour myColour(red, green, blue, alpha);

    G4VModelCommand<M>::Model()->SetDefault(myColour);
  }
}

//Add command
template <typename M>
class G4ModelCommandAdd : public G4VModelCommand<M> {

public: // With description

  G4ModelCommandAdd(M* model, const G4String& placement);

  virtual ~G4ModelCommandAdd();

  void SetNewValue(G4UIcommand* command, G4String newValue);

private:

  G4UIcommand* cmd;
};

template <typename M>
G4ModelCommandAdd<M>::G4ModelCommandAdd(M* model, const G4String& placement)
  :G4VModelCommand<M>(model)
{
  //Set colour through a string
  G4String name = model->Name();
  G4String stringDir = placement+"/"+name+"/add";
  G4UIparameter* param(0);

  cmd = new G4UIcommand(stringDir, this);
  cmd->SetGuidance("Add Command");
  cmd->SetGuidance("One input is expected.");

  param = new G4UIparameter("Value ", 's', false);
  cmd->SetParameter(param);
}

template <typename M>
G4ModelCommandAdd<M>::~G4ModelCommandAdd()
{  
  delete cmd;
}

template <typename M>
void G4ModelCommandAdd<M>::SetNewValue(G4UIcommand* cmd, G4String newValue)
{
  G4VModelCommand<M>::Model()->Add(newValue);
}


//Invert command
template <typename M>
class G4ModelCommandInvert : public G4VModelCommand<M> {

public: // With description

  G4ModelCommandInvert(M* model, const G4String& placement);

  virtual ~G4ModelCommandInvert();

  void SetNewValue(G4UIcommand* command, G4String newValue);

private:

  G4UIcmdWithABool* fpCmd;
};

template <typename M>
G4ModelCommandInvert<M>::G4ModelCommandInvert(M* model, const G4String& placement)
  :G4VModelCommand<M>(model)
{
  //Set colour through a string
  G4String name = model->Name();
  G4String stringDir = placement+"/"+name+"/invert";

  fpCmd = new G4UIcmdWithABool(stringDir, this);
  fpCmd->SetGuidance("Invert command.");
  fpCmd->SetGuidance("One input is expected.");
  fpCmd->SetParameterName("invert", true);
  fpCmd->SetDefaultValue(true);
}

template <typename M>
G4ModelCommandInvert<M>::~G4ModelCommandInvert()
{  
  delete fpCmd;
}

template <typename M>
void G4ModelCommandInvert<M>::SetNewValue(G4UIcommand* cmd, G4String newValue)
{
  G4VModelCommand<M>::Model()->SetInvert(fpCmd->GetNewBoolValue(newValue));
}

//Active command
template <typename M>
class G4ModelCommandActive : public G4VModelCommand<M> {

public: // With description

  G4ModelCommandActive(M* model, const G4String& placement);

  virtual ~G4ModelCommandActive();

  void SetNewValue(G4UIcommand* command, G4String newValue);

private:

  G4UIcmdWithABool* fpCmd;
};

template <typename M>
G4ModelCommandActive<M>::G4ModelCommandActive(M* model, const G4String& placement)
  :G4VModelCommand<M>(model)
{
  //Set colour through a string
  G4String name = model->Name();
  G4String stringDir = placement+"/"+name+"/active";

  fpCmd = new G4UIcmdWithABool(stringDir, this);
  fpCmd->SetGuidance("Active command.");
  fpCmd->SetParameterName("active", true);
  fpCmd->SetDefaultValue(true);
}

template <typename M>
G4ModelCommandActive<M>::~G4ModelCommandActive()
{  
  delete fpCmd;
}

template <typename M>
void G4ModelCommandActive<M>::SetNewValue(G4UIcommand* cmd, G4String newValue)
{
  G4VModelCommand<M>::Model()->SetActive(fpCmd->GetNewBoolValue(newValue));
}


//Reset command
template <typename M>
class G4ModelCommandReset : public G4VModelCommand<M> {

public: // With description

  G4ModelCommandReset(M* model, const G4String& placement);

  virtual ~G4ModelCommandReset();

  void SetNewValue(G4UIcommand* command, G4String newValue);

private:

  G4UIcommand* fpCmd;
};

template <typename M>
G4ModelCommandReset<M>::G4ModelCommandReset(M* model, const G4String& placement)
  :G4VModelCommand<M>(model)
{
  //Set colour through a string
  G4String name = model->Name();
  G4String stringDir = placement+"/"+name+"/reset";

  fpCmd = new G4UIcommand(stringDir, this);
  fpCmd->SetGuidance("Reset command.");
}

template <typename M>
G4ModelCommandReset<M>::~G4ModelCommandReset()
{  
  delete fpCmd;
}

template <typename M>
void G4ModelCommandReset<M>::SetNewValue(G4UIcommand* cmd, G4String string)
{
  G4VModelCommand<M>::Model()->Reset();
}

//Verbose command
template <typename M>
class G4ModelCommandVerbose : public G4VModelCommand<M> {

public: // With description

  G4ModelCommandVerbose(M* model, const G4String& placement);

  virtual ~G4ModelCommandVerbose();

  void SetNewValue(G4UIcommand* command, G4String newValue);

private:

  G4UIcmdWithABool* fpCmd;
};

template <typename M>
G4ModelCommandVerbose<M>::G4ModelCommandVerbose(M* model, const G4String& placement)
  :G4VModelCommand<M>(model)
{
  G4String name = model->Name();
  G4String stringDir = placement+"/"+name+"/verbose";

  fpCmd = new G4UIcmdWithABool(stringDir, this);
  fpCmd->SetGuidance("Verbose command.");
  fpCmd->SetParameterName("verbose", true);
  fpCmd->SetDefaultValue(true);
}

template <typename M>
G4ModelCommandVerbose<M>::~G4ModelCommandVerbose()
{  
  delete fpCmd;
}

template <typename M>
void G4ModelCommandVerbose<M>::SetNewValue(G4UIcommand* cmd, G4String newValue)
{
  G4VModelCommand<M>::Model()->SetVerbose(fpCmd->GetNewBoolValue(newValue));
}

#endif
