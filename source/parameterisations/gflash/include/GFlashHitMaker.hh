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
//
//
//---------------------------------------------------------------
//  GEANT 4 class header file
//
//  GFlashHitMaker
//
//  Class description:
//
//  Singleton for GFlash parameterisation hit creation.

//
// Author: Joanna Weng - 9.11.04
//---------------------------------------------------------------
#ifndef GFlashHitMaker_h
#define GFlashHitMaker_h 1

#include "G4TouchableHandle.hh"
#include "G4Navigator.hh"

#include "GFlashEnergySpot.hh"
#include "G4GFlashSpot.hh"
#include "G4FastTrack.hh"

class G4Step;
class G4StepPoint;
class G4VProcess;

class GFlashHitMaker 
{
  public:

    GFlashHitMaker();
    ~GFlashHitMaker();
  
    void make(GFlashEnergySpot * aSpot, const G4FastTrack * aT );
    inline void SetNameOfWorldWithSD(const G4String& aName) {fWorldWithSdName = aName;};
  
    inline void SetProcess(G4VProcess* proc) { fpProcess = proc; }

  private:  

    G4TouchableHandle fTouchableHandle;
    G4Navigator *fpNavigator;
    G4bool fNaviSetup;
    /// Name of the world containing the sensitive detector. If empty, default mass world is used.
    G4String fWorldWithSdName;

    G4Step* fpSpotS;
    G4StepPoint* fpSpotP;
    G4VProcess* fpProcess = nullptr;

  private:

    GFlashHitMaker(const GFlashHitMaker & ) {}
    GFlashHitMaker & operator = (const GFlashHitMaker & )
    {
      return *this;
    }
};
#endif

