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
// G4GenericMessenger
//
// Author:
//    P.Mato, CERN - 27 September 2012
// Updates:
//    M.Asai, SLAC - 26 November 2013 
//      Adding methods with unit declaration and making thread-safe for
//      version 10.
//    M.Asai, SLAC - 04 May 2014
//      Fix core dump when GetCurrentValue() method is invoked for
//      a command defined by DeclareMethod().
//    M.Asai, SLAC - 11 July 2021
//      Adding G4ThreeVector type without unit
// --------------------------------------------------------------------

#include "G4GenericMessenger.hh"
#include "G4Types.hh"
#include "G4UImessenger.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIdirectory.hh"
#include "G4Threading.hh"

#include <iostream>

class G4InvalidUICommand : public std::bad_cast
{
  public:
   G4InvalidUICommand() = default;
   const char* what() const throw() override
   {
     return "G4InvalidUICommand: command does not exist or is of invalid type";
    }
};

G4GenericMessenger::G4GenericMessenger(void* obj, const G4String& dir,
                                       const G4String& doc)
  : directory(dir)
  , object(obj)
{
  dircmd = new G4UIdirectory(dir);
  dircmd->SetGuidance(doc);
}

G4GenericMessenger::~G4GenericMessenger()
{
  delete dircmd;
  for(const auto& propertie : properties)
  {
    delete propertie.second.command;
  }
  for(const auto& method : methods)
  {
    delete method.second.command;
  }
}

G4GenericMessenger::Command& G4GenericMessenger::DeclareProperty(
  const G4String& name, const G4AnyType& var, const G4String& doc)
{
  G4String fullpath = directory + name;
  G4UIcommand* cmd = nullptr;
  if(var.TypeInfo() == typeid(G4ThreeVector))
  {
    cmd = new G4UIcmdWith3Vector(fullpath.c_str(), this);
    (static_cast<G4UIcmdWith3Vector*>(cmd))
      ->SetParameterName("valueX", "valueY", "valueZ", false, false);
  }
  else
  {
    cmd  = new G4UIcommand(fullpath.c_str(), this);
    char ptype;
    if(var.TypeInfo() == typeid(int) || var.TypeInfo() == typeid(long) ||
       var.TypeInfo() == typeid(unsigned int) ||
       var.TypeInfo() == typeid(unsigned long))
    { ptype = 'i'; }
    else if(var.TypeInfo() == typeid(float) || var.TypeInfo() == typeid(double))
    { ptype = 'd'; }
    else if(var.TypeInfo() == typeid(bool))
    { ptype = 'b'; }
    else if(var.TypeInfo() == typeid(G4String))
    { ptype = 's'; }
    else
    { ptype = 's'; }
    cmd->SetParameter(new G4UIparameter("value", ptype, false));
  }
  if(!doc.empty())
  { cmd->SetGuidance(doc); }
  return properties[name] = Property(var, cmd);
}

G4GenericMessenger::Command& G4GenericMessenger::DeclarePropertyWithUnit(
  const G4String& name, const G4String& defaultUnit, const G4AnyType& var,
  const G4String& doc)
{
  if(var.TypeInfo() != typeid(float) && var.TypeInfo() != typeid(double) &&
     var.TypeInfo() != typeid(G4ThreeVector))
  {
    return DeclareProperty(name, var, doc);
  }
  G4String fullpath = directory + name;
  G4UIcommand* cmd;
  if(var.TypeInfo() == typeid(float) || var.TypeInfo() == typeid(double))
  {
    cmd = new G4UIcmdWithADoubleAndUnit(fullpath.c_str(), this);
    (static_cast<G4UIcmdWithADoubleAndUnit*>(cmd))
      ->SetParameterName("value", false, false);
    (static_cast<G4UIcmdWithADoubleAndUnit*>(cmd))->SetDefaultUnit(defaultUnit);
  }
  else
  {
    cmd = new G4UIcmdWith3VectorAndUnit(fullpath.c_str(), this);
    (static_cast<G4UIcmdWith3VectorAndUnit*>(cmd))
      ->SetParameterName("valueX", "valueY", "valueZ", false, false);
    (static_cast<G4UIcmdWith3VectorAndUnit*>(cmd))->SetDefaultUnit(defaultUnit);
  }

  if(!doc.empty())
  {
    cmd->SetGuidance(doc);
  }
  return properties[name] = Property(var, cmd);
}

G4GenericMessenger::Command& G4GenericMessenger::DeclareMethod(
  const G4String& name, const G4AnyMethod& fun, const G4String& doc)
{
  G4String fullpath = directory + name;
  auto* cmd         = new G4UIcommand(fullpath.c_str(), this);
  if(!doc.empty())
  {
    cmd->SetGuidance(doc);
  }
  for(std::size_t i = 0; i < fun.NArg(); ++i)
  {
    G4String argNam = "arg" + ItoS((G4int)i);
    char ptype = 's';
    auto& tInfo = fun.ArgType(i);
    if(tInfo == typeid(int) || tInfo == typeid(long) ||
       tInfo == typeid(unsigned int) ||
       tInfo == typeid(unsigned long))
    { ptype = 'i'; }
    else if(tInfo == typeid(float) || tInfo == typeid(double))
    { ptype = 'd'; }
    else if(tInfo == typeid(bool))
    { ptype = 'b'; }
    else if(tInfo == typeid(G4String))
    { ptype = 's'; }
    else
    { ptype = 's'; }
    cmd->SetParameter(new G4UIparameter(argNam, ptype, false));
  }
  return methods[name] = Method(fun, object, cmd);
}

G4GenericMessenger::Command& G4GenericMessenger::DeclareMethodWithUnit(
  const G4String& name, const G4String& defaultUnit, const G4AnyMethod& fun,
  const G4String& doc)
{
  G4String fullpath = directory + name;
  if(fun.NArg() != 1)
  {
    G4ExceptionDescription ed;
    ed << "G4GenericMessenger::DeclareMethodWithUnit() does not support a "
          "method that has more than\n"
       << "one arguments (or no argument). Please use "
          "G4GenericMessenger::DeclareMethod method for\n"
       << "your command <" << fullpath << ">.";
    G4Exception("G4GenericMessenger::DeclareMethodWithUnit()", "Intercom70002",
                FatalException, ed);
  }
  G4UIcommand* cmd = new G4UIcmdWithADoubleAndUnit(fullpath.c_str(), this);
  (static_cast<G4UIcmdWithADoubleAndUnit*>(cmd))
    ->SetParameterName("value", false, false);
  (static_cast<G4UIcmdWithADoubleAndUnit*>(cmd))->SetDefaultUnit(defaultUnit);
  if(!doc.empty())
  {
    cmd->SetGuidance(doc);
  }
  return methods[name] = Method(fun, object, cmd);
}

G4String G4GenericMessenger::GetCurrentValue(G4UIcommand* command)
{
  if(properties.find(command->GetCommandName()) != properties.cend())
  {
    Property& p = properties[command->GetCommandName()];
    return p.variable.ToString();
  }
  else if(methods.find(command->GetCommandName()) != methods.cend())
  {
    G4cout << " GetCurrentValue() is not available for a command defined by "
              "G4GenericMessenger::DeclareMethod()."
           << G4endl;
    return G4String();
  }
  else
  {
    throw G4InvalidUICommand();
  }
}

void G4GenericMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  // Check if there are units on this commands
  if(typeid(*command) == typeid(G4UIcmdWithADoubleAndUnit))
  {
    newValue = G4UIcommand::ConvertToString(
      G4UIcommand::ConvertToDimensionedDouble(newValue));
  }
  else if(typeid(*command) == typeid(G4UIcmdWith3VectorAndUnit))
  {
    newValue = G4UIcommand::ConvertToString(
      G4UIcommand::ConvertToDimensioned3Vector(newValue));
  }

  if(properties.find(command->GetCommandName()) != properties.cend())
  {
    Property& p = properties[command->GetCommandName()];
    p.variable.FromString(newValue);
  }
  else if(methods.find(command->GetCommandName()) != methods.cend())
  {
    Method& m = methods[command->GetCommandName()];
    if(m.method.NArg() == 0)
    {
      m.method.operator()(m.object);
    }
    else if(m.method.NArg() > 0)
    {
      m.method.operator()(m.object, newValue);
    }
    else
    {
      throw G4InvalidUICommand();
    }
  }
}

void G4GenericMessenger::SetGuidance(const G4String& s)
{
  dircmd->SetGuidance(s);
}

G4GenericMessenger::Command& G4GenericMessenger::Command::SetUnit(
  const G4String& unit, UnitSpec spec)
{
  // Change the type of command (unfortunatelly this is done a posteriory)
  // We need to delete the old command before creating the new one and therefore
  // we need to recover the information before the deletetion
  if(G4Threading::IsMultithreadedApplication())
  {
    G4String cmdpath = command->GetCommandPath();
    G4ExceptionDescription ed;
    ed << "G4GenericMessenger::Command::SetUnit() is thread-unsafe and should "
          "not be used\n"
       << "in multi-threaded mode. For your command <" << cmdpath << ">, use\n"
       << " DeclarePropertyWithUnit(const G4String& name, const G4String& "
          "defaultUnit,\n"
       << "                         const G4AnyType& variable, const G4String& "
          "doc)\n"
       << "or\n"
       << " DeclareMethodWithUnit(const G4String& name, const G4String& "
          "defaultUnit,\n"
       << "                       const G4AnyType& variable, const G4String& "
          "doc)\n"
       << "to define a command with a unit <" << unit << ">.";
    if(spec != UnitDefault)
    {
      ed << "\nPlease use a default unit instead of unit category.";
    }
    G4Exception("G4GenericMessenger::Command::SetUnit()", "Intercom70001",
                FatalException, ed);
    return *this;
  }

  G4String cmdpath         = command->GetCommandPath();
  G4UImessenger* messenger = command->GetMessenger();
  G4String range           = command->GetRange();
  std::vector<G4String> guidance;
  G4String par_name   = command->GetParameter(0)->GetParameterName();
  G4bool par_omitable = command->GetParameter(0)->IsOmittable();
  for(G4int i = 0; i < (G4int)command->GetGuidanceEntries(); ++i)
  {
    guidance.push_back(command->GetGuidanceLine(i));
  }
  // Before deleting the command we need to add a fake one to avoid deleting
  // the directory entry and with its guidance
  G4UIcommand tmp((cmdpath + "_tmp").c_str(), messenger);
  delete command;

  if(*type == typeid(float) || *type == typeid(double))
  {
    auto* cmd_t = new G4UIcmdWithADoubleAndUnit(cmdpath, messenger);
    if(spec == UnitDefault)
    {
      cmd_t->SetDefaultUnit(unit);
    }
    else if(spec == UnitCategory)
    {
      cmd_t->SetUnitCategory(unit);
    }
    cmd_t->SetParameterName(par_name, par_omitable);
    command = cmd_t;
  }
  else if(*type == typeid(G4ThreeVector))
  {
    auto* cmd_t = new G4UIcmdWith3VectorAndUnit(cmdpath, messenger);
    if(spec == UnitDefault)
    {
      cmd_t->SetDefaultUnit(unit);
    }
    else if(spec == UnitCategory)
    {
      cmd_t->SetUnitCategory(unit);
    }
    command = cmd_t;
  }
  else
  {
    G4cerr << "Only parameters of type <double> or <float> can be associated "
              "with units"
           << G4endl;
    return *this;
  }
  for(auto& i : guidance)
  {
    command->SetGuidance(i);
  }
  command->SetRange(range);
  return *this;
}

G4GenericMessenger::Command& G4GenericMessenger::Command::SetParameterName(
  const G4String& name, G4bool omittable, G4bool currentAsDefault)
{ 
  return SetParameterName(0,name,omittable,currentAsDefault);
}

G4GenericMessenger::Command& G4GenericMessenger::Command::SetParameterName(
  G4int pIdx,
  const G4String& name, G4bool omittable, G4bool currentAsDefault)
{
  if(pIdx < 0 || pIdx >= (G4int)(command->GetParameterEntries()))
  {
    G4cerr << "Invalid parameter index : " << pIdx << "\nMethod ignored."
           << G4endl;
    return *this;
  }
  G4UIparameter* theParam = command->GetParameter(pIdx);
  theParam->SetParameterName(name);
  theParam->SetOmittable(omittable);
  theParam->SetCurrentAsDefault(currentAsDefault);
  return *this;
}

G4GenericMessenger::Command& G4GenericMessenger::Command::SetParameterName(
  const G4String& namex, const G4String& namey, const G4String& namez,
  G4bool omittable, G4bool currentAsDefault)
{
  if(*type != typeid(G4ThreeVector))
  {
    G4cerr << "This SetParameterName method is for G4ThreeVector!! " 
           << "Method ignored." << G4endl;
    return *this;
  }
  G4UIparameter* theParam = command->GetParameter(0);
  theParam->SetParameterName(namex);
  theParam->SetOmittable(omittable);
  theParam->SetCurrentAsDefault(currentAsDefault);
  theParam = command->GetParameter(1);
  theParam->SetParameterName(namey);
  theParam->SetOmittable(omittable);
  theParam->SetCurrentAsDefault(currentAsDefault);
  theParam = command->GetParameter(2);
  theParam->SetParameterName(namez);
  theParam->SetOmittable(omittable);
  theParam->SetCurrentAsDefault(currentAsDefault);
  return *this;
}

G4GenericMessenger::Command& G4GenericMessenger::Command::SetCandidates(
  const G4String& candList)
{
  return SetCandidates(0,candList);
}

G4GenericMessenger::Command& G4GenericMessenger::Command::SetCandidates(
  G4int pIdx, const G4String& candList)
{
  if(pIdx < 0 || pIdx >= (G4int)(command->GetParameterEntries()))
  {
    G4cerr << "Invalid parameter index : " << pIdx << "\nMethod ignored."
           << G4endl;
    return *this;
  }
  G4UIparameter* theParam = command->GetParameter(pIdx);
  theParam->SetParameterCandidates(candList);
  return *this;
}

G4GenericMessenger::Command& G4GenericMessenger::Command::SetDefaultValue(
  const G4String& defVal)
{
  return SetDefaultValue(0,defVal);
}

G4GenericMessenger::Command& G4GenericMessenger::Command::SetDefaultValue(
  G4int pIdx, const G4String& defVal)
{
  if(pIdx < 0 || pIdx >= (G4int)(command->GetParameterEntries()))
  {
    G4cerr << "Invalid parameter index : " << pIdx << "\nMethod ignored."
           << G4endl;
    return *this;
  }
  G4UIparameter* theParam = command->GetParameter(pIdx);
  theParam->SetDefaultValue(defVal);
  return *this;
}
