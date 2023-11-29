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

#ifndef G4FASTSIMHITMAKER_HH
#define G4FASTSIMHITMAKER_HH

#include "G4TouchableHandle.hh"
#include "G4Navigator.hh"
#include "G4FastHit.hh"
#include "G4FastTrack.hh"
class G4Step;
class G4StepPoint;
class G4VProcess;

/**
 * @brief Helper class for hit creation
 *
 * Helper class that can be employed in the fast simulation models.
 * It allows to deposit energy at given position (G4FastHit), provided it is
 * located within the sensitive detector that derives from 
 * G4VFastSimSensitiveDetector base class.
 * An extended example extended/parameterisations/Par03 demonstrates how to use
 * G4FastSimHitMaker to create multiple deposits from the fast simulation model.
 *
 */

class G4FastSimHitMaker
{
 public:
  G4FastSimHitMaker();
  ~G4FastSimHitMaker();

  /// Deposit energy at given position.
  /// @param[in] aHit Created hit (energy and position)
  /// @param[in] aTrack Fast track with access to particle's track and
  /// properties in envelope's local coordinates
  void make(const G4FastHit& aHit, const G4FastTrack& aTrack);
  /// If sensitive detector class is in the parallel world, it must be
  /// specified, otherwise no sensitive detector will be found (mass geometry
  /// will be checked).
  /// @param[in] aName Name of the parallel world
  inline void SetNameOfWorldWithSD(const G4String& aName)
  {
    fWorldWithSdName = aName;
  };
  inline void SetProcess(G4VProcess* proc) { fpProcess = proc; }

 private:
  /// Touchable
  G4TouchableHandle fTouchableHandle;
  /// Navigator
  G4Navigator* fpNavigator;
  /// Flag specifying if navigator has been already set up
  G4bool fNaviSetup;
  /// Name of the world containing the sensitive detector. If empty, default
  /// mass world is used.
  G4String fWorldWithSdName;

  G4Step* fpSpotS;
  G4StepPoint* fpSpotP;
  G4VProcess* fpProcess = nullptr;
};
#endif
