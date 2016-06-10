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
// $Id: G4UImessenger.hh,v 1.9 2006-06-29 19:08:19 gunter Exp $
//

#ifndef G4GenericMessenger_h
#define G4GenericMessenger_h 1

#include "G4UImessenger.hh"
#include "G4UIcommand.hh"
#include "G4AnyType.hh"
#include "G4AnyMethod.hh"
#include "G4ApplicationState.hh"

#include <map>
#include <vector>

class G4UIdirectory;

/// This class is generic messenger.

class G4GenericMessenger : public G4UImessenger
{
public:
  /// Contructor
  G4GenericMessenger(void* obj, const G4String& dir = "", const G4String& doc = "");
  /// Destructor
  virtual ~G4GenericMessenger();
  /// The concrete, but generic implementation of this method.
  virtual G4String GetCurrentValue(G4UIcommand* command);
  /// The concrete, generic  implementation of this method converts the string "newValue" to action.
  virtual void SetNewValue(G4UIcommand* command, G4String newValue);
  
public:
  struct Command {
    enum UnitSpec {UnitCategory, UnitDefault};
    Command(G4UIcommand* cmd, const std::type_info& ti) : command(cmd), type(&ti) {}
    Command() : command(0), type(0) {}
    //    Command& operator =(const Command& rhs) { command = rhs.command; type = rhs.type; }
    Command& SetStates(G4ApplicationState s0) {command->AvailableForStates(s0); return *this;}
    Command& SetStates(G4ApplicationState s0, G4ApplicationState s1) {command->AvailableForStates(s0, s1); return *this;}
    Command& SetStates(G4ApplicationState s0, G4ApplicationState s1, G4ApplicationState s2){command->AvailableForStates(s0,s1,s2); return *this;}
    Command& SetStates(G4ApplicationState s0, G4ApplicationState s1, G4ApplicationState s2, G4ApplicationState s3) {command->AvailableForStates(s0,s1,s2,s3); return *this;}
    Command& SetStates(G4ApplicationState s0, G4ApplicationState s1, G4ApplicationState s2, G4ApplicationState s3, G4ApplicationState s4) {command->AvailableForStates(s0,s1,s2,s3,s4); return *this;}
    Command& SetRange(const G4String& range) {command->SetRange(range.c_str()); return *this;}
    Command& SetGuidance(const G4String& s0) { command->SetGuidance(s0); return *this; }
    Command& SetUnit(const G4String&, UnitSpec = UnitDefault);
    Command& SetUnitCategory(const G4String& u) {return SetUnit(u, UnitCategory);}
    Command& SetDefaultUnit(const G4String& u) {return SetUnit(u, UnitDefault);}
    Command& SetParameterName(const G4String&, G4bool, G4bool =false);
    Command& SetDefaultValue(const G4String&);
    Command& SetCandidates(const G4String&);
    Command& SetToBeBroadcasted(G4bool s0) { command->SetToBeBroadcasted(s0); return *this; }
    Command& SetToBeFlushed(G4bool s0) { command->SetToBeFlushed(s0); return *this; }
    Command& SetWorkerThreadOnly(G4bool s0) { command->SetWorkerThreadOnly(s0); return *this; }
    
    G4UIcommand* command;
    const std::type_info* type;
  };
  struct Property : public Command {
    Property(const G4AnyType& var, G4UIcommand* cmd) : Command(cmd, var.TypeInfo()) , variable(var) {}
    Property() {}
    G4AnyType variable;
  };
  struct Method : public Command {
    Method(const G4AnyMethod& fun, void* obj, G4UIcommand* cmd) : Command(cmd, fun.ArgType()), method(fun), object(obj) {}
    Method() : object(0) {}
    G4AnyMethod method;
    void* object;
  };
  
  ///Declare Methods
  Command& DeclareProperty(const G4String& name, const G4AnyType& variable, const G4String& doc = "");
  Command& DeclarePropertyWithUnit
   (const G4String& name, const G4String& defaultUnit, const G4AnyType& variable, const G4String& doc = "");
  Command& DeclareMethod(const G4String& name, const G4AnyMethod& fun, const G4String& doc = "");
  Command& DeclareMethodWithUnit
   (const G4String& name, const G4String& defaultUnit, const G4AnyMethod& fun, const G4String& doc = "");
  void SetDirectory(const G4String& dir) {directory = dir;}
  void SetGuidance(const G4String& s);
    
private:
  std::map<G4String, Property> properties;
  std::map<G4String, Method> methods;
  G4UIdirectory* dircmd;
  G4String directory;
  void* object;
};


#endif

