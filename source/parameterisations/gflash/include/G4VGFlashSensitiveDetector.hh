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
// $Id: G4VGFlashSensitiveDetector.hh 97674 2016-06-07 08:28:08Z gcosmo $
//
//
//---------------------------------------------------------------
//  GEANT 4 class header file
//
//  G4VGFlashSensitiveDetector
//
//  Class description:
//
// Abstract base class of the sensitive detector for use with GFlash.
// The user's sensitive detector which generates hits must be derived
// from this class, and G4VSensitiveDetector.

//---------------------------------------------------------------
#ifndef G4VGFlashSensitiveDetector_h
#define G4VGFlashSensitiveDetector_h 1

#include "G4Step.hh"
#include "G4VReadOutGeometry.hh"
#include "G4TouchableHistory.hh"
#include "GFlashEnergySpot.hh"
#include "G4GFlashSpot.hh"
#include "G4VSensitiveDetector.hh"


class G4VGFlashSensitiveDetector 
{

  public: // with description

      G4VGFlashSensitiveDetector() {}
      G4VGFlashSensitiveDetector(const G4VGFlashSensitiveDetector &) {}
       // Constructors. The user's concrete class must use one of these
       // constructors by the constructor initializer of the derived class.
       // The name of the sensitive detector must be the same as for the
       // corresponding GG4VSensitiveDetector.

  public: // without description

      virtual ~G4VGFlashSensitiveDetector() {}

      G4int operator==(const G4VGFlashSensitiveDetector &right) const
        {return this == &right;}
      G4int operator!=(const G4VGFlashSensitiveDetector &right) const
        {return this != &right;}

  public: // without description

      inline G4bool Hit(G4GFlashSpot * aSpot)
      {
        // This is the public method invoked by GFlashHitMaker for generating
        // hits. The actual user's implementation for generating hits must be
        // implemented in GenerateHits() virtual protected method. 

        G4bool result = true; 
        G4VSensitiveDetector * This
          = dynamic_cast<G4VSensitiveDetector *>(this);
        if(!This)
        {
          G4Exception("G4VGFlashSensitiveDetector::Hit()",
                      "InvalidSetup", FatalException,
                      "Needs also to inherit from G4VSensitiveDetector!");
          return false;
        }
        if(This->isActive())
        { 
          G4VReadOutGeometry * ROgeometry = 0;
          G4TouchableHistory* ROhis = 0;

          if(This) ROgeometry = This->GetROgeometry();
          if(ROgeometry)
          {
            // fake pre-step point for touchable from read-out geometry.
            G4Step fakeStep;
            G4StepPoint * tmpPoint = fakeStep.GetPreStepPoint();
            tmpPoint->SetTouchableHandle(aSpot->GetTouchableHandle());
            tmpPoint->SetPosition(aSpot->GetPosition());
            tmpPoint->SetMomentumDirection(aSpot->GetOriginatorTrack()
                               ->GetPrimaryTrack()->GetMomentumDirection());
            result = ROgeometry->CheckROVolume(&fakeStep, ROhis); 
          }
          if(result) result = ProcessHits(aSpot, ROhis); 
        }
        else 
        {
          result = false;
        }
        return result;
      }

  protected: // with description

      virtual G4bool ProcessHits(G4GFlashSpot*aSpot,
                                 G4TouchableHistory*ROhist) = 0;
       // The user MUST implement this method for generating hit(s) from the
       // GFlashSpots. Be aware that this method is a protected method and it
       // will be invoked by Hit() method of the Base class once the Readout
       // geometry that may be associated to the corresponding
       // G4VSensitiveDetector was taken into account.
};

#endif

