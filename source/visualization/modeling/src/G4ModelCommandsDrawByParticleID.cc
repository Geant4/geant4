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
// $Id: G4ModelCommandsDrawByParticleID.cc,v 1.3 2005/11/28 20:07:11 tinslay Exp $
// GEANT4 tag $Name: geant4-08-00 $
// 
// Jane Tinslay, John Allison, Joseph Perl November 2005

#include "G4ModelCommandsDrawByParticleID.hh"
#include "G4TrajectoryDrawByParticleID.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcommand.hh"
#include <sstream>

//Set colour with string
G4ModelCommandDrawByParticleIDSet::G4ModelCommandDrawByParticleIDSet(G4TrajectoryDrawByParticleID* model, const G4String& placement)
  :G4VModelCommand<G4TrajectoryDrawByParticleID>(model) 
{
  G4String name = model->Name();
  G4String myCommand = placement+"/"+name+"/set";

  fpCommand = new G4UIcmdWithAString(myCommand, this);      
  fpCommand->SetGuidance("Set trajectory colour through a string.");
  fpCommand->SetGuidance("Two inputs are expected, for example");

  G4String example = myCommand+" gamma red";

  fpCommand->SetGuidance(example);
  fpCommand->SetParameterName("parameters", false); 
}

G4ModelCommandDrawByParticleIDSet::~G4ModelCommandDrawByParticleIDSet() 
{
  delete fpCommand;
} 

void G4ModelCommandDrawByParticleIDSet::SetNewValue(G4UIcommand*, G4String newValue) 
{
  G4String particle;
  G4String colour;
  std::istringstream is (newValue);
  is >> particle >> colour;  

  // Configure model as requested
  Model()->Set(particle, colour);
}

//Set colour by red, green, blue and alpha components
G4ModelCommandDrawByParticleIDSetRGBA::G4ModelCommandDrawByParticleIDSetRGBA(G4TrajectoryDrawByParticleID* model, const G4String& placement) 
  :G4VModelCommand<G4TrajectoryDrawByParticleID>(model)
{
  G4String name = model->Name();
  G4String myCommand = placement+"/"+name+"/setRGBA";

  fpCommand = new G4UIcmdWithAString(myCommand, this);      
  fpCommand->SetGuidance("Set trajectory colour through red, green, blue and alpha components.");
  fpCommand->SetGuidance("Five inputs are expected, for example");

  G4String example = myCommand+" gamma 1 1 1 1";
;
  fpCommand->SetGuidance(example);
  fpCommand->SetParameterName("parameters", false); 

}

G4ModelCommandDrawByParticleIDSetRGBA::~G4ModelCommandDrawByParticleIDSetRGBA() 
{
  delete fpCommand;
} 

void G4ModelCommandDrawByParticleIDSetRGBA::SetNewValue(G4UIcommand*, G4String newValue) 
{
  G4String particle;
  G4double red(0), green(0), blue(0), alpha(0);
  std::istringstream is (newValue);
  is >> particle >> red >> green >> blue >> alpha;  

  G4Colour myColour(red, green, blue, alpha);
  
  // Configure model as requested
  Model()->Set(particle, myColour);
}


//Set default colour with string
G4ModelCommandDrawByParticleIDSetDefault::G4ModelCommandDrawByParticleIDSetDefault(G4TrajectoryDrawByParticleID* model, const G4String& placement)
  :G4VModelCommand<G4TrajectoryDrawByParticleID>(model) 
{
  G4String name = model->Name();
  G4String myCommand = placement+"/"+name+"/setDefault";

  fpCommand = new G4UIcmdWithAString(myCommand, this);      
  fpCommand->SetGuidance("Set default trajectory colour through a string.");
  fpCommand->SetGuidance("One input is expected, for example");

  G4String example = myCommand+" red";

  fpCommand->SetGuidance(example);
  fpCommand->SetParameterName("parameters", false); 
}

G4ModelCommandDrawByParticleIDSetDefault::~G4ModelCommandDrawByParticleIDSetDefault() 
{
  delete fpCommand;
} 

void G4ModelCommandDrawByParticleIDSetDefault::SetNewValue(G4UIcommand*, G4String newValue) 
{
  // Configure model as requested
  Model()->SetDefault(newValue);
}

//Set colour by red, green, blue and alpha components
G4ModelCommandDrawByParticleIDSetDefaultRGBA::G4ModelCommandDrawByParticleIDSetDefaultRGBA(G4TrajectoryDrawByParticleID* model, const G4String& placement) 
  :G4VModelCommand<G4TrajectoryDrawByParticleID>(model)
{
  G4String name = model->Name();
  G4String myCommand = placement+"/"+name+"/setDefaultRGBA";

  fpCommand = new G4UIcmdWithAString(myCommand, this);      
  fpCommand->SetGuidance("Set default trajectory colour through red, green, blue and alpha components.");
  fpCommand->SetGuidance("Four inputs are expected, for example");

  G4String example = myCommand+" 1 1 1 1";
;
  fpCommand->SetGuidance(example);
  fpCommand->SetParameterName("parameters", false); 
}

G4ModelCommandDrawByParticleIDSetDefaultRGBA::~G4ModelCommandDrawByParticleIDSetDefaultRGBA() 
{
  delete fpCommand;
} 

void G4ModelCommandDrawByParticleIDSetDefaultRGBA::SetNewValue(G4UIcommand*, G4String newValue) 
{
  G4double red(0), green(0), blue(0), alpha(0);
  std::istringstream is (newValue);
  is >> red >> green >> blue >> alpha;  

  G4Colour myColour(red, green, blue, alpha);
  
  // Configure model as requested
  Model()->SetDefault(myColour);
}
