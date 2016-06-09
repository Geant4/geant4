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
// File name:     RadmonMaterialsMessenger.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonMaterialsMessenger.cc,v 1.4 2006/06/29 16:16:53 gunter Exp $
// Tag:           $Name: geant4-08-01 $
//

// Messenger commands path
#define COMMANDS_PATH "/radmon/materials/"

// Include files
#include "RadmonMaterialsMessenger.hh"
#include "RadmonMaterialsManager.hh"
#include "G4UnitsTable.hh"



                                                RadmonMaterialsMessenger :: RadmonMaterialsMessenger(RadmonMaterialsManager * manager)
:
 RadmonMessenger(COMMANDS_PATH, "Interactive materials definition commands."),
 materialsManager(manager),
 RADMON_INITIALIZE_COMMAND(CreateElement),
 RADMON_INITIALIZE_COMMAND(CreateMaterial),
 RADMON_INITIALIZE_COMMAND(AddComponentByAtoms),
 RADMON_INITIALIZE_COMMAND(AddComponentByFraction),
 RADMON_INITIALIZE_COMMAND(SetMaterialColor),
 RADMON_INITIALIZE_COMMAND(SetMaterialTrasparency),
 RADMON_INITIALIZE_COMMAND(SetMaterialVisibility),
 RADMON_INITIALIZE_COMMAND(SetMaterialStyle),
 RADMON_INITIALIZE_COMMAND(Dump),
 RADMON_INITIALIZE_COMMAND(Insert),
 RADMON_INITIALIZE_COMMAND(Save)
{
 RADMON_CREATE_COMMAND_5ARGS(CreateElement,             "Creates a new element",                                                        "material", "symbol", "Z", "A", "unit");
 RADMON_CREATE_COMMAND_4ARGS(CreateMaterial,            "Creates a new material",                                                       "material", "density", "unit", "totComponents");
 RADMON_CREATE_COMMAND_3ARGS(AddComponentByAtoms,       "Adds a element into a material specifying the number of atoms",                "material", "element", "atoms");
 RADMON_CREATE_COMMAND_3ARGS(AddComponentByFraction,    "Adds a element or material into another material specifying the percentage",   "material", "component", "fraction");
 RADMON_CREATE_COMMAND_4ARGS(SetMaterialColor,          "Set default color for all elements made of this material",                     "material", "red", "green", "blue");
 RADMON_CREATE_COMMAND_2ARGS(SetMaterialTrasparency,    "Set default trasparency for all elements made of this material",               "material", "alpha");
 RADMON_CREATE_COMMAND_2ARGS(SetMaterialVisibility,     "Set default visibility for all elements made of this material",                "material", "visibility");
 RADMON_CREATE_COMMAND_2ARGS(SetMaterialStyle,          "Set default style for all elements made of this material",                     "material", "style");
 RADMON_CREATE_COMMAND_0ARGS(Dump,                      "Print currently declared elements and materials");
 RADMON_CREATE_COMMAND_1ARG (Insert,                    "Inserts elements and materials from file",                                     "fileName");
 RADMON_CREATE_COMMAND_1ARG (Save,                      "Saves a elements and materials declarations into a file",                      "fileName");
}



                                                RadmonMaterialsMessenger :: ~RadmonMaterialsMessenger()
{
 RADMON_DESTROY_COMMAND(Save);
 RADMON_DESTROY_COMMAND(Insert);
 RADMON_DESTROY_COMMAND(Dump);
 RADMON_DESTROY_COMMAND(SetMaterialStyle);
 RADMON_DESTROY_COMMAND(SetMaterialVisibility);
 RADMON_DESTROY_COMMAND(SetMaterialTrasparency);
 RADMON_DESTROY_COMMAND(SetMaterialColor);
 RADMON_DESTROY_COMMAND(AddComponentByFraction);
 RADMON_DESTROY_COMMAND(AddComponentByAtoms);
 RADMON_DESTROY_COMMAND(CreateMaterial);
 RADMON_DESTROY_COMMAND(CreateElement);
}





G4String                                        RadmonMaterialsMessenger :: GetCurrentValue(G4UIcommand * /* command */)
{
 G4cout << "RadmonMaterialsMessenger::GetCurrentValue(): Not supported" << G4endl;
 
 return G4String();
}



void                                            RadmonMaterialsMessenger :: SetNewValue(G4UIcommand * command, G4String newValue)
{
 RADMON_BEGIN_LIST_SET_COMMANDS
  RADMON_SET_COMMAND(CreateElement)
  RADMON_SET_COMMAND(CreateMaterial)
  RADMON_SET_COMMAND(AddComponentByAtoms)
  RADMON_SET_COMMAND(AddComponentByFraction)
  RADMON_SET_COMMAND(SetMaterialColor)
  RADMON_SET_COMMAND(SetMaterialTrasparency)
  RADMON_SET_COMMAND(SetMaterialVisibility)
  RADMON_SET_COMMAND(SetMaterialStyle)
  RADMON_SET_COMMAND(Dump)
  RADMON_SET_COMMAND(Insert)
  RADMON_SET_COMMAND(Save)
 RADMON_END_LIST_SET_COMMANDS
}





// Events
void                                            RadmonMaterialsMessenger :: OnCreateElement(const G4String & value)
{
 G4String args[5];

 if (!ProcessArguments(value, 5, args))
  return; 

 G4double z(G4UIcommand::ConvertToDouble(args[2]));
 
 if (z<=0.)
 {
  G4cout << "RadmonMaterialsMessenger::OnCreateElement(): z must be positive." << G4endl;
  return;
 }
 
 G4double a(GetUnit(args[4], "Molar Mass"));
 
 if (a<0.)
  return;
 
 a*=G4UIcommand::ConvertToDouble(args[3]);
 
 if (a<=0.)
 {
  G4cout << "RadmonMaterialsMessenger::OnCreateElement(): a must be positive." << G4endl;
  return;
 }

 materialsManager->CreateElement(args[0], args[1], z, a);
}



void                                            RadmonMaterialsMessenger :: OnCreateMaterial(const G4String & value)
{
 G4String args[4];

 if (!ProcessArguments(value, 4, args))
  return; 

 G4double density(GetUnit(args[2], "Volumic Mass"));
 
 if (density<0.)
  return;
 
 density*=G4UIcommand::ConvertToDouble(args[1]);
 
 if (density<=0.)
 {
  G4cout << "RadmonMaterialsMessenger::OnCreateMaterial(): Material density must be positive." << G4endl;
  return;
 }
 
 G4int components(G4UIcommand::ConvertToInt(args[3]));
 
 if (components<0)
 {
  G4cout << "RadmonMaterialsMessenger::OnCreateMaterial(): Number of components must be positive (>=0)." << G4endl;
  return;
 }

 materialsManager->CreateMaterial(args[0], density, components);
}



void                                            RadmonMaterialsMessenger :: OnAddComponentByAtoms(const G4String & value)
{
 G4String args[3];

 if (!ProcessArguments(value, 3, args))
  return; 

 if (!materialsManager->IsIncompleteMaterial(args[0]))
 {
  if (!materialsManager->ExistsMaterial(args[0]))
   G4cout << "RadmonMaterialsMessenger::OnAddComponentByAtoms(): Material \"" << args[0] << "\" not found." << G4endl;
  else
   G4cout << "RadmonMaterialsMessenger::OnAddComponentByAtoms(): No more components to be added to \"" << args[0] << "\"." << G4endl;
  return;
 }

 if (!materialsManager->ExistsElement(args[1]))
 {
  G4cout << "RadmonMaterialsMessenger::OnAddComponentByAtoms(): Element \"" << args[1] << "\" not found." << G4endl;
  return;
 }

 G4int atoms(G4UIcommand::ConvertToInt(args[2]));
 
 if (atoms<=0)
 {
  G4cout << "RadmonMaterialsMessenger::OnAddComponentByAtoms(): Number of atoms must be positive (>0)." << G4endl;
  return;
 }

 materialsManager->AddComponentByAtoms(args[0], args[1], atoms); 
}



void                                            RadmonMaterialsMessenger :: OnAddComponentByFraction(const G4String & value)
{
 G4String args[3];

 if (!ProcessArguments(value, 3, args))
  return; 

 if (!materialsManager->IsIncompleteMaterial(args[0]))
 {
  if (!materialsManager->ExistsMaterial(args[0]))
   G4cout << "RadmonMaterialsMessenger::OnAddComponentByFraction(): Material \"" << args[0] << "\" not found." << G4endl;
  else
   G4cout << "RadmonMaterialsMessenger::OnAddComponentByFraction(): No more components to be added to \"" << args[0] << "\"." << G4endl;
  return;
 }

 if ((!materialsManager->ExistsElement(args[1])) && (!materialsManager->ExistsMaterial(args[1])))
 {
  G4cout << "RadmonMaterialsMessenger::OnAddComponentByFraction(): Element/Material \"" << args[1] << "\" not found." << G4endl;
  return;
 }

 G4double fraction(G4UIcommand::ConvertToDouble(args[2]));
 
 if (fraction<=0)
 {
  G4cout << "RadmonMaterialsMessenger::OnAddComponentByFraction(): Material fraction must be positive (>0)." << G4endl;
  return;
 }

 materialsManager->AddComponentByFraction(args[0], args[1], fraction); 
}



void                                            RadmonMaterialsMessenger :: OnSetMaterialColor(const G4String & value)
{
 G4String args[4];
 
 if (!ProcessArguments(value, 4, args))
  return;
  
 if (!materialsManager->ExistsMaterial(args[0]) && !materialsManager->IsIncompleteMaterial(args[0]))
 {
  G4cout << "RadmonMaterialsMessenger::OnSetMaterialColor(): Material \"" << args[0] << "\" not found." << G4endl;
  return;
 }
 
 G4double red(G4UIcommand::ConvertToDouble(args[1]));
 G4double green(G4UIcommand::ConvertToDouble(args[2]));
 G4double blue(G4UIcommand::ConvertToDouble(args[3]));
 
 if (red<0. || red>1. || green<0. || green>1. || blue<0. || blue>1.)
 {
  G4cout << "RadmonMaterialsMessenger::OnSetMaterialColor(): Red, green and blue components must be set between 0 and 1." << G4endl;
  return;
 }
 
 materialsManager->SetMaterialColor(args[0], G4Color(red, green, blue, materialsManager->GetMaterialColor(args[0]).GetAlpha()));
}



void                                            RadmonMaterialsMessenger :: OnSetMaterialTrasparency(const G4String & value)
{
 G4String args[2];
 
 if (!ProcessArguments(value, 2, args))
  return;
  
 if (!materialsManager->ExistsMaterial(args[0]) && !materialsManager->IsIncompleteMaterial(args[0]))
 {
  G4cout << "RadmonMaterialsMessenger::OnSetMaterialTrasparency(): Material \"" << args[0] << "\" not found." << G4endl;
  return;
 }
 
 G4double alpha(G4UIcommand::ConvertToDouble(args[1]));
 
 if (alpha<0. || alpha>1.)
 {
  G4cout << "RadmonMaterialsMessenger::OnSetMaterialTrasparency(): Alpha component must be set between 0 and 1." << G4endl;
  return;
 }
 
 const G4Color &c(materialsManager->GetMaterialColor(args[0]));
 materialsManager->SetMaterialColor(args[0], G4Color(c.GetRed(), c.GetGreen(), c.GetBlue(), alpha));
}



void                                            RadmonMaterialsMessenger :: OnSetMaterialVisibility(const G4String & value)
{
 G4String args[2];
 
 if (!ProcessArguments(value, 2, args))
  return;
  
 if (!materialsManager->ExistsMaterial(args[0]) && !materialsManager->IsIncompleteMaterial(args[0]))
 {
  G4cout << "RadmonMaterialsMessenger::OnSetMaterialVisibility(): Material \"" << args[0] << "\" not found." << G4endl;
  return;
 }
 
 if (args[1]=="hidden")
  materialsManager->SetMaterialVisibility(args[0], false);
 else if (args[1]=="visible")
  materialsManager->SetMaterialVisibility(args[0], true);
 else
  G4cout << "RadmonMaterialsMessenger::OnSetMaterialVisibility(): Visibility must be set to visible or hidden." << G4endl;
}



void                                            RadmonMaterialsMessenger :: OnSetMaterialStyle(const G4String & value)
{
 G4String args[2];
 
 if (!ProcessArguments(value, 2, args))
  return;
  
 if (!materialsManager->ExistsMaterial(args[0]) && !materialsManager->IsIncompleteMaterial(args[0]))
 {
  G4cout << "RadmonMaterialsMessenger::OnSetMaterialStyle(): Material \"" << args[0] << "\" not found." << G4endl;
  return;
 }
 
 if (args[1]=="wireframe")
  materialsManager->SetMaterialForceWireframe(args[0], true);
 else if (args[1]=="solid")
  materialsManager->SetMaterialForceSolid(args[0], true);
 else if (args[1]=="default")
 {
  materialsManager->SetMaterialForceWireframe(args[0], false);
  materialsManager->SetMaterialForceSolid(args[0], false);
 }
 else
  G4cout << "RadmonMaterialsMessenger::OnSetMaterialStyle(): Style must be set to wireframe, solid or default." << G4endl;
}



void                                            RadmonMaterialsMessenger :: OnDump(const G4String & /* value */)
{
 materialsManager->Dump(G4cout, "- ");
 G4cout << G4endl;
}



void                                            RadmonMaterialsMessenger :: OnInsert(const G4String & value)
{
 G4String args[1];

 if (!ProcessArguments(value, 1, args))
  return; 
  
 std::istream * in(OpenForInput(args[0]));
 
 if (!in)
  return;

 if (!materialsManager->Insert(*in)) 
  G4cout << "RadmonMaterialsMessenger::OnInsert(): Error reading from file \"" << args[0] << "\"." << G4endl;
  
 delete in;
}



void                                            RadmonMaterialsMessenger :: OnSave(const G4String & value)
{
 G4String args[1];

 if (!ProcessArguments(value, 1, args))
  return; 
  
 std::ostream * out(OpenForOutput(args[0]));
 
 if (!out)
  return;

 if (!materialsManager->Save(*out))
  G4cout << "RadmonMaterialsMessenger::OnSave(): Cannot write layout into file \"" << args[0] << "\"." << G4endl;
  
 delete out;
}

