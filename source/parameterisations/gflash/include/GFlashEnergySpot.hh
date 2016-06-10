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
// $Id: GFlashEnergySpot.hh 68057 2013-03-13 14:46:00Z gcosmo $
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
