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

#ifndef G4VFASTSIMSENSITIVEDETECTOR_HH
#define G4VFASTSIMSENSITIVEDETECTOR_HH

#include "G4VReadOutGeometry.hh"
#include "G4TouchableHistory.hh"
#include "G4VSensitiveDetector.hh"
#include "G4FastHit.hh"
#include "G4FastTrack.hh"

/**
 * @brief Base class for the sensitive detector used within the fast simulation
 *
 * Base class for a sensitive detector that allows to store hits created in the
 * fast simulation models. It must me used in addition to inheritance from the
 * usual base class `G4VSensitiveDetector` for processing of energy deposited
 * in G4Step.
 * ProcessHits(...) method must be implemented and describe how hits should be
 * saved in the hit collections. It is invoked by Hit method which is public and
 * can be called directly in the fast simulation model (if sensitive detector is
 * known to the model), or via the helper class G4FastSimHitMaker that will
 * locate appropriate volume and retrieve its sensitive detector.
 * An extended example extended/parameterisations/Par03 demonstrates how to use
 * G4VFastSimSensitiveDetector to deposit energy from fast simulation and
 * compare it to the detailed simulation.
 *
 */

class G4VFastSimSensitiveDetector
{
 public:
  /// Create a hit.
  ///
  /// It checks if G4VSensitiveDetector is also used as a base class,
  /// and takes into account the readout geometry, if it is defined.
  /// User instruction on how to deposit energy needs to be implemented in
  /// ProcessHits method.
  /// @param[in] aHit Created hit (energy and position)
  /// @param[in] aTrack Fast track with access to particle's track and
  /// properties in envelope's local coordinates
  /// @param[in] aTouchable Touchable with relevant transformations
  inline G4bool Hit(const G4FastHit* aHit, const G4FastTrack* aTrack,
                    G4TouchableHandle* aTouchable)
  {
    G4bool result                 = true;
    G4VSensitiveDetector* sensDet = dynamic_cast<G4VSensitiveDetector*>(this);
    if(sensDet == nullptr)
    {
      G4Exception("G4VFastSimSensitiveDetector::Hit()", "InvalidSetup",
                  FatalException,
                  "Sensitive detector needs also to inherit also from "
                  "G4VSensitiveDetector if full "
                  "simulation is used instead!");
    }
    if(sensDet->isActive())
    {
      G4VReadOutGeometry* ROgeometry = sensDet->GetROgeometry();
      G4TouchableHistory* ROhistory  = 0;

      if(ROgeometry)
      {
        // create fake pre-step point updating the touchable from read-out
        // geometry.
        G4Step fakeStep;
        const G4Track* currentTrack = aTrack->GetPrimaryTrack();
        G4StepPoint* tmpPoint       = fakeStep.GetPreStepPoint();
        tmpPoint->SetTouchableHandle(*aTouchable);
        tmpPoint->SetPosition(aHit->GetPosition());
        tmpPoint->SetMomentumDirection(currentTrack->GetMomentumDirection());
        result = ROgeometry->CheckROVolume(&fakeStep, ROhistory);
      } else {
        ROhistory = static_cast<G4TouchableHistory*>((*aTouchable)());
      }
      if(result)
        result = ProcessHits(aHit, aTrack, ROhistory);
    }
    else
    {
      result = false;
    }
    return result;
  }

 private:
  /// Describes how energy and position of deposits are inserted into the hits
  /// collection. It is a private method and it will be invoked by Hit() method
  /// of the base class once the readout geometry that may be associated to the
  /// corresponding G4VSensitiveDetector is taken into account.
  /// It needs to be implemented in the derived class.
  /// @param[in] aHit Created hit (energy and position)
  /// @param[in] aTrack Fast track with access to particle's track and
  /// properties in envelope's local coordinates
  /// @param[in] aROHistory Touchable history with relevant transformations
  virtual G4bool ProcessHits(const G4FastHit* aHit, const G4FastTrack* aTrack,
                             G4TouchableHistory* aROHistory) = 0;
};
#endif /* G4VFASTSIMSENSITIVEDETECTOR_HH */
