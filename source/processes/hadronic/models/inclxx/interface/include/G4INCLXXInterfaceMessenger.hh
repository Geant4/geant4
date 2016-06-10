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
// Alain Boudard, CEA-Saclay, France
// Joseph Cugnon, University of Liege, Belgium
// Jean-Christophe David, CEA-Saclay, France
// Pekka Kaitaniemi, CEA-Saclay, France, and Helsinki Institute of Physics, Finland
// Sylvie Leray, CEA-Saclay, France
// Davide Mancusi, CEA-Saclay, France
//
#define INCLXX_IN_GEANT4_MODE 1

#include "globals.hh"

/** \file G4INCLXXInterfaceMessenger.hh
 * \brief Messenger class for the Geant4 INCL++ interface.
 *
 * \date 26th April 2012
 * \author Davide Mancusi
 */

#ifndef G4INCLXXInterfaceMessenger_hh
#define G4INCLXXInterfaceMessenger_hh

#include "G4INCLXXInterfaceStore.hh"
#include "G4UImessenger.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4String.hh"

class G4INCLXXInterfaceStore;

class G4INCLXXInterfaceMessenger : public G4UImessenger
{

  public:
    G4INCLXXInterfaceMessenger (G4INCLXXInterfaceStore *anInterfaceStore);
    ~G4INCLXXInterfaceMessenger ();
    void SetNewValue (G4UIcommand *command, G4String newValues);
    //    Identifies the command which has been invoked by the user, extracts the
    //    parameters associated with that command (held in newValues, and uses
    //    these values with the appropriate member function of G4INCLXXInterfaceStore.
    //
  private:
    static const G4String theUIDirectory;
    G4INCLXXInterfaceStore *theINCLXXInterfaceStore;
    G4UIdirectory *theINCLXXDirectory;
    G4UIcmdWithAString *accurateNucleusCmd;
    G4UIcmdWithAnInteger *maxClusterMassCmd;
    G4UIcmdWithADoubleAndUnit *cascadeMinEnergyPerNucleonCmd;
    G4UIcmdWithAString *inclPhysicsCmd;
    G4UIcommand *useAblaCmd;
};

#endif

