//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: GFlashHitMaker.hh,v 1.4 2005/10/04 09:08:33 gcosmo Exp $
// GEANT4 tag $Name: geant4-08-00 $
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

class GFlashHitMaker 
{
  public:

    GFlashHitMaker();
    ~GFlashHitMaker();
  
    void make(GFlashEnergySpot * aSpot, const G4FastTrack * aT );
  
  private:  

    G4TouchableHandle fTouchableHandle;
    G4Navigator *fpNavigator;
    G4bool fNaviSetup;
  
  private:

    GFlashHitMaker(const GFlashHitMaker & ) {}
    GFlashHitMaker & operator = (const GFlashHitMaker & )
    {
      return *this;
    }
};
#endif

