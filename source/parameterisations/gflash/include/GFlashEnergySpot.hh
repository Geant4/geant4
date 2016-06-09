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
// $Id: GFlashEnergySpot.hh,v 1.4 2005/10/04 09:08:33 gcosmo Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
//
//---------------------------------------------------------------
//  GEANT 4 class header file
//
//  GFlashEnergySpot
//
//  Class description:
//
//  Alternative definition of energy spot for GFlash parameterisation.

//
// Author: Joanna Weng - 9.11.04
//---------------------------------------------------------------
#ifndef GFlashEnergySpot_h
#define GFlashEnergySpot_h

#include "G4ThreeVector.hh"

class GFlashEnergySpot
{
  public:

    GFlashEnergySpot();
    GFlashEnergySpot(const G4ThreeVector& point, G4double E);
    ~GFlashEnergySpot();
  
    inline void SetEnergy(const G4double& E) {Energy = E;}
    inline G4double GetEnergy() const {return Energy;}
  
    inline void SetPosition(const G4ThreeVector& point) {Point = point;}
    inline G4ThreeVector GetPosition() const {return Point;}
    
  private:

    G4double Energy;     // energy deposition
    G4ThreeVector Point; // locus of energy deposition
};

#endif
