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
// INCL++ intra-nuclear cascade model
// Pekka Kaitaniemi, CEA and Helsinki Institute of Physics
// Davide Mancusi, CEA
// Alain Boudard, CEA
// Sylvie Leray, CEA
// Joseph Cugnon, University of Liege
//
// INCL++ revision: v5.1.4
//
#define INCLXX_IN_GEANT4_MODE 1

#include "globals.hh"

/** \file G4INCLXXInterfaceConfig.hh
 * \brief Singleton class for configuring the INCL++ Geant4 interface.
 *
 * \date 24 May 2012
 * \author Davide Mancusi
 */

#ifndef G4INCLXXINTERFACECONFIG_HH_
#define G4INCLXXINTERFACECONFIG_HH_

#include "G4INCLXXInterface.hh"
#include <list>

class G4INCLXXInterfaceMessenger;

class G4INCLXXInterfaceConfig {
  // Typedefs to simplify some code
  typedef std::list<G4INCLXXInterface *> G4INCLXXInterfaceList;
  typedef std::list<G4INCLXXInterface *>::const_iterator G4INCLXXInterfaceIter;

  public:

  /// \brief Get the singleton instance
  static G4INCLXXInterfaceConfig *GetInstance() {
    if(!theInstance)
      return (theInstance = new G4INCLXXInterfaceConfig);
    else
      return theInstance;
  }

  /// \brief Delete the singleton instance
  static void DeleteInstance() {
    delete theInstance;
    theInstance = NULL;
  }

  void RegisterINCLXXInterface(G4INCLXXInterface * const anInterface) {
    theInterfaces.push_back(anInterface);
  }

  void DeleteModels() {
    for(G4INCLXXInterfaceIter anInterface=theInterfaces.begin();
        anInterface!=theInterfaces.end();
        ++anInterface) {
      (*anInterface)->DeleteModel();
    }
  }

  /// \brief Setter for accurateProjectile
  void SetAccurateProjectile(const G4bool b) { accurateProjectile=b; }

  /// \brief Setter for theMaxClusterMass
  void SetMaxClusterMass(const G4int aMass) { theMaxClusterMass=aMass; }




  /// \brief Getter for accurateProjectile
  G4bool GetAccurateProjectile() { return accurateProjectile; }

  /// \brief Getter for ClusterMaxMass
  G4int GetMaxClusterMass() { return theMaxClusterMass; }

  /// \brief Getter for dumpInput
  G4bool GetDumpInput() { return dumpInput; }

  /// \brief Getter for theMaxProjMassINCL
  G4int GetMaxProjMassINCL() { return theMaxProjMassINCL; }

  private:
  /** \brief Private constructor
   *
   * This class is a singleton. It must be instantiated using the GetInstance
   * static method.
   */
  G4INCLXXInterfaceConfig();

  /** \brief Private destructor
   *
   * This class is a singleton. Its instance must be deleted using the
   * DeleteInstance static method.
   */
  ~G4INCLXXInterfaceConfig();

  static G4INCLXXInterfaceConfig *theInstance;

  G4bool dumpInput;
  G4bool accurateProjectile;
  const G4int theMaxClusterMassDefault;
  G4int theMaxClusterMass;
  const G4int theMaxProjMassINCL;

  G4INCLXXInterfaceMessenger *theINCLXXInterfaceMessenger;

  G4INCLXXInterfaceList theInterfaces;
};

#endif // G4INCLXXINTERFACECONFIG_HH_
