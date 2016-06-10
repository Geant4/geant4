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
// $Id: G4GFlashSpot.hh 68057 2013-03-13 14:46:00Z gcosmo $
//
//
//---------------------------------------------------------------
//  GEANT 4 class header file
//
//  G4GFlashSpot
//
//  Class description:
//
//  Definition of energy spot for GFlash parameterisation.

//---------------------------------------------------------------
#ifndef G4GFlashSpot_h
#define G4GFlashSpot_h

#include "G4ThreeVector.hh"
#include "G4FastTrack.hh"
#include "GFlashEnergySpot.hh"
#include "G4TouchableHandle.hh"

class G4GFlashSpot
{
  public:

    G4GFlashSpot(const GFlashEnergySpot * aSpot,
                 const G4FastTrack * aTrack, G4TouchableHandle aH)
      : theSpot(aSpot), theTrack(aTrack), theHandle(aH) {}
  
    ~G4GFlashSpot() {}
  
    const GFlashEnergySpot * GetEnergySpot() const {return theSpot;}
  
    const G4FastTrack * GetOriginatorTrack() const {return theTrack;}
  
    G4TouchableHandle GetTouchableHandle() const {return theHandle;}
  
    G4ThreeVector GetPosition() const
     {return GetOriginatorTrack()->GetPrimaryTrack()->GetPosition();}
    
  private:
  
    const GFlashEnergySpot * theSpot;
    const G4FastTrack * theTrack;
    G4TouchableHandle theHandle;
};

#endif
