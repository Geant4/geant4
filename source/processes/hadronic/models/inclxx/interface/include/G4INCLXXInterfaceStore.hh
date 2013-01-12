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
#define INCLXX_IN_GEANT4_MODE 1

#include "globals.hh"

/** \file G4INCLXXInterfaceStore.hh
 * \brief Header file for the G4INCLXXInterfaceStore class
 *
 * \date 24 May 2012
 * \author Davide Mancusi
 */

#ifndef G4INCLXXINTERFACESTORE_HH_
#define G4INCLXXINTERFACESTORE_HH_

#include "G4INCLXXInterface.hh"
#include "G4INCLCascade.hh"
#include "G4INCLConfig.hh"
#include <list>
#include <sstream>

class G4INCLXXInterfaceMessenger;

/** \class G4INCLXXInterfaceStore
 * \brief Singleton class for configuring the INCL++ Geant4 interface.
 *
 * This class also contains a single cached instance of the INCL model
 * (\see{G4INCL::INCL}).
 */
class G4INCLXXInterfaceStore {
  public:

    /// \brief Get the singleton instance
    static G4INCLXXInterfaceStore *GetInstance() {
      if(!theInstance)
        theInstance = new G4INCLXXInterfaceStore;
      return theInstance;
    }

    /// \brief Delete the singleton instance
    static void DeleteInstance() {
      delete theInstance;
      theInstance = NULL;
    }

    /// \brief Get the cached INCL model engine
    G4INCL::INCL *GetINCLModel() {
      if(!theINCLModel) {
        G4INCL::Config *theConfig = new G4INCL::Config;
        theConfig->setClusterMaxMass(theMaxClusterMass);
        theINCLModel = new G4INCL::INCL(theConfig);
        // ownership of the Config object is taken over by the INCL model engine
      }
      return theINCLModel;
    }




    /// \brief Setter for accurateProjectile
    void SetAccurateProjectile(const G4bool b) {
      if(accurateProjectile!=b) {
        // Parameter is changed, emit a big warning message
        std::stringstream ss;
        ss << "Switching from "
          << (accurateProjectile ? "\"accurate projectile\" mode to \"accurate target\"" : "\"accurate target\" mode to \"accurate projectile\"")
          << " mode."
          << G4endl
          << "Do this ONLY if you fully understand what it does!";
        EmitBigWarning(ss.str());
      }

      // No need to delete the model for this parameter

      accurateProjectile=b;
    }

    /// \brief Setter for theMaxClusterMass
    void SetMaxClusterMass(const G4int aMass) {
      if(theMaxClusterMass!=aMass) {
        // Parameter is changed, emit a big warning message
        std::stringstream ss;
        ss << "Changing maximum cluster mass from "
          << theMaxClusterMass
          << " to "
          << aMass
          << "."
          << G4endl
          << "Do this ONLY if you fully understand what this setting does!";
        EmitBigWarning(ss.str());
      }

      // We must delete the model object to make sure that we use the new
      // parameter
      DeleteModel();

      theMaxClusterMass=aMass;
    }




    /** \brief Getter for accurateProjectile
     *
     * The \see{G4INCLXXInterfaceMessenger} class provides a UI command to set
     * this parameter.
     */
    G4bool GetAccurateProjectile() const { return accurateProjectile; }

    /** \brief Getter for ClusterMaxMass
     *
     * The \see{G4INCLXXInterfaceMessenger} class provides a UI command to set
     * this parameter.
     */
    G4int GetMaxClusterMass() const { return theMaxClusterMass; }




    /// \brief Getter for theMaxProjMassINCL
    G4int GetMaxProjMassINCL() const { return theMaxProjMassINCL; }

    /// \brief Getter for dumpInput
    G4bool GetDumpInput() const { return dumpInput; }





    /** \brief Emit a warning to G4cout
     *
     * The InterfaceStore will not emit more than maxWarnings warnings.
     */
    void EmitWarning(const G4String &message);

    /** \brief Emit a BIG warning to G4cout
     *
     * There is no limit on the number of BIG warnings emitted.
     */
    void EmitBigWarning(const G4String &message) const;

  private:

    /** \brief Private constructor
     *
     * This class is a singleton. It must be instantiated using the GetInstance
     * static method.
     */
    G4INCLXXInterfaceStore();

    /** \brief Private destructor
     *
     * This class is a singleton. Its instance must be deleted using the
     * DeleteInstance static method.
     */
    ~G4INCLXXInterfaceStore();

    /// \brief Delete the INCL model engine
    void DeleteModel() { delete theINCLModel; theINCLModel=NULL; }

    /// \brief Create a new Config object from the current options

    static G4INCLXXInterfaceStore *theInstance;

    G4bool dumpInput;
    G4bool accurateProjectile;
    const G4int theMaxClusterMassDefault;
    G4int theMaxClusterMass;
    const G4int theMaxProjMassINCL;

    G4INCLXXInterfaceMessenger *theINCLXXInterfaceMessenger;

    G4INCL::INCL *theINCLModel;

    /// \brief Static warning counter
    G4int nWarnings;

    /// \brief Maximum number of warnings
    const G4int maxWarnings;
};

#endif // G4INCLXXINTERFACESTORE_HH_
