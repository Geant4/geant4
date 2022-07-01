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
// G4GenericMessenger
//
// Class description:
//
// A generic messenger class.

// Author: P.Mato, CERN - 27 September 2012
// --------------------------------------------------------------------
#ifndef G4GenericMessenger_hh
#define G4GenericMessenger_hh 1

#include "G4UImessenger.hh"
#include "G4UIcommand.hh"
#include "G4AnyType.hh"
#include "G4AnyMethod.hh"
#include "G4ApplicationState.hh"

#include <map>
#include <vector>

class G4UIdirectory;

class G4GenericMessenger : public G4UImessenger
{
  public:

    G4GenericMessenger(void* obj, const G4String& dir = "",
                       const G4String& doc = "");
      // Contructor

    ~G4GenericMessenger() override;
    // Destructor

    G4String GetCurrentValue(G4UIcommand* command) override;
    // The concrete, but generic implementation of this method.

    void SetNewValue(G4UIcommand* command, G4String newValue) override;
    // The concrete, generic  implementation of this method converts
    // the string "newValue" to action.

   public:

    struct Command
    {
      enum UnitSpec
      {
        UnitCategory,
        UnitDefault
      };
      Command(G4UIcommand* cmd, const std::type_info& ti)
        : command(cmd)
        , type(&ti)
      {}
      Command() = default;

      Command& SetStates(G4ApplicationState s0)
      {
        command->AvailableForStates(s0);
        return *this;
      }
      Command& SetStates(G4ApplicationState s0, G4ApplicationState s1)
      {
        command->AvailableForStates(s0, s1);
        return *this;
      }
      Command& SetStates(G4ApplicationState s0, G4ApplicationState s1,
                         G4ApplicationState s2)
      {
        command->AvailableForStates(s0, s1, s2);
        return *this;
      }
      Command& SetStates(G4ApplicationState s0, G4ApplicationState s1,
                         G4ApplicationState s2, G4ApplicationState s3)
      {
        command->AvailableForStates(s0, s1, s2, s3);
        return *this;
      }
      Command& SetStates(G4ApplicationState s0, G4ApplicationState s1,
                         G4ApplicationState s2, G4ApplicationState s3,
                         G4ApplicationState s4)
      {
        command->AvailableForStates(s0, s1, s2, s3, s4);
        return *this;
      }
      Command& SetRange(const G4String& range)
      {
        command->SetRange(range.c_str());
        return *this;
      }
      Command& SetGuidance(const G4String& s0)
      {
        command->SetGuidance(s0);
        return *this;
      }
      Command& SetUnit(const G4String&, UnitSpec = UnitDefault);
      Command& SetUnitCategory(const G4String& u)
      {
        return SetUnit(u, UnitCategory);
      }
      Command& SetDefaultUnit(const G4String& u)
      {
        return SetUnit(u, UnitDefault);
      }
      Command& SetParameterName(const G4String&, G4bool, G4bool = false);
      Command& SetParameterName(G4int pIdx, const G4String&, G4bool, G4bool = false);
      Command& SetParameterName(const G4String&, const G4String&, const G4String&, 
                                G4bool, G4bool = false);
      Command& SetDefaultValue(const G4String&);
      Command& SetDefaultValue(G4int pIdx, const G4String&);
      Command& SetCandidates(const G4String&);
      Command& SetCandidates(G4int pIdx, const G4String&);
      Command& SetToBeBroadcasted(G4bool s0)
      {
        command->SetToBeBroadcasted(s0);
        return *this;
      }
      Command& SetToBeFlushed(G4bool s0)
      {
        command->SetToBeFlushed(s0);
        return *this;
      }
      Command& SetWorkerThreadOnly(G4bool s0)
      {
        command->SetWorkerThreadOnly(s0);
        return *this;
      }

      G4UIcommand* command = nullptr;
      const std::type_info* type = nullptr;
    };

    struct Property : public Command
    {
      Property(const G4AnyType& var, G4UIcommand* cmd)
        : Command(cmd, var.TypeInfo())
        , variable(var)
      {}
      Property() = default;
      G4AnyType variable;
    };

    struct Method : public Command
    {
      Method(const G4AnyMethod& fun, void* obj, G4UIcommand* cmd)
        : Command(cmd, fun.ArgType())
        , method(fun)
        , object(obj)
      {}
      Method() = default;
      G4AnyMethod method;
      void* object = nullptr;
    };

    // Declare Methods

    Command& DeclareProperty(const G4String& name, const G4AnyType& variable,
                             const G4String& doc = "");
    Command& DeclarePropertyWithUnit(const G4String& name,
                                     const G4String& defaultUnit,
                                     const G4AnyType& variable,
                                     const G4String& doc = "");
    Command& DeclareMethod(const G4String& name, const G4AnyMethod& fun,
                           const G4String& doc = "");
    Command& DeclareMethodWithUnit(const G4String& name,
                                   const G4String& defaultUnit,
                                   const G4AnyMethod& fun,
                                   const G4String& doc = "");
    void SetDirectory(const G4String& dir) { directory = dir; }
    void SetGuidance(const G4String& s);
    void Sort(G4bool val = true)
    {
      if(dircmd != nullptr)
      {
        dircmd->Sort(val);
      }
    }

  private:

    std::map<G4String, Property> properties;
    std::map<G4String, Method> methods;
    G4UIdirectory* dircmd = nullptr;
    G4String directory;
    void* object = nullptr;
};

#endif
