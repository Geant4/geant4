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
// $Id:
//
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:    G4PhysListRegistry
//
// Author  R. Hatcher  2014-10-15
//
// Modifications:  based on G4PhysicsConstructorRegistry
//
//
// Class Description
// This is a singleton keeping pointers to all physics lists
// Class Description - End

#ifndef G4PhysListRegistry_hh
#define G4PhysListRegistry_hh 1

#include <vector>
#include <map>
#include "globals.hh"


class G4VModularPhysicsList;
class G4VBasePhysListStamper;

class G4PhysListRegistry
{
public:

  static G4PhysListRegistry* Instance();
  // access

  ~G4PhysListRegistry();

  void AddFactory(G4String name, G4VBasePhysListStamper*);
  void AddPhysicsExtension(G4String name, G4String procname);
  // mapping from extension name to actual physics contructor process name

  G4VModularPhysicsList* GetModularPhysicsList(const G4String& name);
  G4VModularPhysicsList* GetModularPhysicsListFromEnv();

  G4bool IsReferencePhysList(G4String nam) const;

  const std::vector<G4String>&  AvailablePhysLists() const;
  const std::vector<G4String>&  AvailablePhysicsExtensions() const;
  const std::vector<G4String>&  AvailablePhysListsEM() const;

  void PrintAvailablePhysLists() const;

  G4bool DeconstructPhysListName(const G4String& name, G4String& plBase,
                                 std::vector<G4String>& physExt,
                                 std::vector<G4int>& replace,
                                 G4int verbose=0) const;

  G4bool FindLongestMatch(const G4String& workName,
                          const G4String& searchName,
                          const std::vector<G4String>& validNames,
                          G4String& bestMatch,
                          G4int verbose=0) const;

  inline void  SetVerbose(G4int val) { verbose = val; }
  inline G4int GetVerbose() const { return verbose; }

  inline void  SetUnknownFatal(G4int val) { unknownFatal = val; }
  inline G4int GetUnknownFatal() const { return unknownFatal; }

         void      SetUserDefaultPhysList(const G4String& name="");
  inline G4String  GetUserDefaultPhysList() const { return userDefault; }
  // set a prefered list in case where $PHYSLIST isn't defined
  // if not set (or called with "") this falls back to system default

  inline G4String  GetSystemDefaultPhysList() const { return systemDefault; }

private:

  G4PhysListRegistry();

  static G4ThreadLocal G4PhysListRegistry* theInstance;

  std::map <G4String, G4VBasePhysListStamper*> factories;
  std::map <G4String, G4String> physicsExtensions;

  G4int    verbose;
  G4int    unknownFatal;  /// throw an exception if unsatisfiable?
  G4String userDefault;   /// use this if $PHYSLIST isn't set
  G4String systemDefault; /// use this if user hasn't set userDefault
                          /// or attempts to set the userDefault=""

  // Make these mutable and update them on each request because the map might
  // have been updated by the addition of new entries
  // The only reason to have them at all is that that original interface
  // returned a const reference and so we can't pass back a local object
  // created upon the call of the method
  mutable std::vector<G4String> availBasePhysLists;
  mutable std::vector<G4String> availExtensions;

};

#endif
