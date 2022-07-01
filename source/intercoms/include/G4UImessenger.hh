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
// G4UImessenger
//
// Class description:
//
// This class is the base class representing a messenger which keeps all basic
// commands. The user who wants to define some commands must create his/her
// own concrete class derived from this class. The user's concrete messenger
// must have a responsibility of creating and deleting commands. Also, it must
// take care of the delivering of the commands to the destination class and
// provide the current value(s) of the parameter(s)

// Author: Makoto Asai, 1998
// --------------------------------------------------------------------
#ifndef G4UImessenger_hh
#define G4UImessenger_hh 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4UIdirectory.hh"

class G4UImessenger
{
  public:

    G4UImessenger() = default;
    G4UImessenger(const G4String& path, const G4String& dsc,
                  G4bool commandsToBeBroadcasted = true);
      // Constructor. In the implementation of the concrete messenger,
      // all commands related to the messenger must be constructed

    virtual ~G4UImessenger();
      // Destructor. In the implementation of the concrete messenger,
      // all commands defined in the constructor must be deleted

    virtual G4String GetCurrentValue(G4UIcommand* command);
      // The concrete implementation of this method gets the current value(s)
      // of the parameter(s) of the given command from the destination class,
      // converts the value(s) to a string, and returns the string.
      // Conversion could be done by the ConvertToString() method of
      // corresponding G4UIcmdXXX classes if the command is an object of
      // these G4UIcmdXXX classes

    virtual void SetNewValue(G4UIcommand* command, G4String newValue);
      // The concrete implementation of this method converts the string
      // "newValue" to value(s) of type(s) of the parameter(s).
      // Converted methods corresponding to the type of the command can be
      // used if the command is an object of G4UIcmdXXX classes

    G4bool operator==(const G4UImessenger& messenger) const;
    G4bool operator!=(const G4UImessenger& messenger) const;

    inline G4bool CommandsShouldBeInMaster() const
    {
      return commandsShouldBeInMaster;
    }

  protected:

    G4String ItoS(G4int i);
    G4String DtoS(G4double a);
    G4String BtoS(G4bool b);
    G4int StoI(const G4String& s);
    G4long StoL(const G4String& s);
    G4double StoD(const G4String& s);
    G4bool StoB(G4String s);

    void AddUIcommand(G4UIcommand* newCommand);

    void CreateDirectory(const G4String& path, const G4String& dsc,
                         G4bool commandsToBeBroadcasted = true);
    template <typename T>
    T* CreateCommand(const G4String& cname, const G4String& dsc);
      // Shortcut way for creating directory and commands

  protected:

    G4UIdirectory* baseDir = nullptr;  // used if new object is created
    G4String baseDirName = "";    // used if dir already exists
    G4bool commandsShouldBeInMaster = false;
};

// Inline template implementations

template <typename T>
T* G4UImessenger::CreateCommand(const G4String& cname, const G4String& dsc)
{
  G4String path;
  if(cname[0] != '/')
  {
    path = baseDirName + cname;
    if(path[0] != '/')
    {
      path = "/" + path;
    }
  }

  T* command = new T(path.c_str(), this);
  command->SetGuidance(dsc.c_str());

  return command;
}

#endif
