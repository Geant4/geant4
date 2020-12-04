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
// Author: Jonathan Madsen (May 28th 2020)
//
// class description:
//   This is a class creates a G4RunManager instance. Once can specify an
//   enumerated value to set the default. If "G4RunManagerType::Default"
//   is used,
//

#ifndef G4RunManagerFactory_hh
#define G4RunManagerFactory_hh 1

#include "G4Types.hh"
#include "G4RunManager.hh"
#include "G4MTRunManager.hh"
#include "G4TaskRunManager.hh"
#include "G4VUserTaskQueue.hh"

#include <set>
#include <map>
#include <string>
#include <regex>

//============================================================================//

enum class G4RunManagerType : G4int
{
  Serial      = 0,
  SerialOnly  = 1,
  MT          = 2,
  MTOnly      = 3,
  Tasking     = 4,
  TaskingOnly = 5,
  TBB         = 6,
  TBBOnly     = 7,
  Default
};

//============================================================================//

class G4RunManagerFactory
{
 public:
  static G4RunManager* CreateRunManager(
    G4RunManagerType _type   = G4RunManagerType::Default,
    G4VUserTaskQueue* _queue = nullptr, G4bool fail_if_unavail = true,
    G4int nthreads = 0);

  // provide a version which specifies to fail if unavailable
  static G4RunManager* CreateRunManager(G4RunManagerType _type,
                                        G4bool fail_if_unavail,
                                        G4int nthreads           = 0,
                                        G4VUserTaskQueue* _queue = nullptr)
  {
    return CreateRunManager(_type, _queue, fail_if_unavail, nthreads);
  }

  static G4RunManager* CreateRunManager(G4RunManagerType _type, G4int nthreads,
                                        G4bool fail_if_unavail   = true,
                                        G4VUserTaskQueue* _queue = nullptr)
  {
    return CreateRunManager(_type, _queue, fail_if_unavail, nthreads);
  }

  // provide string versions of above
  template <typename... Args>
  static G4RunManager* CreateRunManager(std::string type, Args&&... args)
  {
    return CreateRunManager(GetType(type), std::forward<Args>(args)...);
  }

  static std::string GetDefault();
  static std::string GetName(G4RunManagerType);
  static G4RunManagerType GetType(const std::string&);
  static std::set<std::string> GetOptions();
  static G4RunManager* GetMasterRunManager();
  static G4MTRunManager* GetMTMasterRunManager();
  static G4RunManagerKernel* GetMasterRunManagerKernel();
};

#endif  // G4TaskRunManagerCreator_h
